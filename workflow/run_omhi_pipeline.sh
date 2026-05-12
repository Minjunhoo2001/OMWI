#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# OMHI full inference pipeline
# QC -> host removal -> MetaPhlAn4 -> merge -> species table -> CLR -> OMHI
# =========================================================

# -----------------------------
# 1. Usage
# -----------------------------
if [ "$#" -lt 2 ]; then
  echo "Usage:"
  echo "  bash workflow/run_omhi_pipeline.sh samples.tsv results_dir"
  echo ""
  echo "samples.tsv format:"
  echo -e "sample_id\tr1\tr2"
  echo -e "sample01\t/path/sample01_R1.fastq.gz\t/path/sample01_R2.fastq.gz"
  echo -e "sample02\t/path/sample02_R1.fastq.gz\t/path/sample02_R2.fastq.gz"
  exit 1
fi

SAMPLE_SHEET="$1"
OUTDIR="$2"

# -----------------------------
# 2. User settings
# -----------------------------
THREADS="${THREADS:-8}"

# Human Bowtie2 index prefix.
# Example:
#   /path/to/GRCh38
# The actual files should look like:
#   /path/to/GRCh38.1.bt2
#   /path/to/GRCh38.2.bt2
#   ...
HOST_INDEX="${HOST_INDEX:-/path/to/GRCh38_bowtie2_index/GRCh38}"

# MetaPhlAn4 database directory and index name
MPA_DB="${MPA_DB:-/share/home/HeMinjun/db/metaphlan4}"
MPA_INDEX="${MPA_INDEX:-mpa_vJun23_CHOCOPhlAnSGB_202403}"

# OMHI model files
OMHI_COEF="${OMHI_COEF:-model/omhi_species_coefficients.tsv}"
OMHI_UNIVERSE="${OMHI_UNIVERSE:-model/omhi_species_universe.tsv}"

# Scripts
LEVEL_EXTRACTOR_SCRIPT="${LEVEL_EXTRACTOR_SCRIPT:-scripts/metaphlan_level_extractor.py}"
CALC_OMHI_SCRIPT="${CALC_OMHI_SCRIPT:-scripts/calc_omhi.R}"

# -----------------------------
# 3. Output folders
# -----------------------------
QC_DIR="${OUTDIR}/01.qc"
NOHOST_DIR="${OUTDIR}/02.nohost"
MPA_DIR="${OUTDIR}/03.metaphlan"
MERGED_DIR="${OUTDIR}/04.merged"
OMHI_DIR="${OUTDIR}/05.omhi"
LOG_DIR="${OUTDIR}/logs"

mkdir -p \
  "${QC_DIR}" \
  "${NOHOST_DIR}" \
  "${MPA_DIR}" \
  "${MERGED_DIR}" \
  "${OMHI_DIR}" \
  "${LOG_DIR}"

MAIN_LOG="${LOG_DIR}/run_omhi_pipeline.log"

exec > >(tee -a "${MAIN_LOG}") 2>&1

echo "[OMHI] Pipeline started at: $(date)"
echo "[OMHI] Sample sheet: ${SAMPLE_SHEET}"
echo "[OMHI] Output directory: ${OUTDIR}"
echo "[OMHI] Threads: ${THREADS}"
echo "[OMHI] Host Bowtie2 index: ${HOST_INDEX}"
echo "[OMHI] MetaPhlAn4 DB: ${MPA_DB}"
echo "[OMHI] MetaPhlAn4 index: ${MPA_INDEX}"
echo "[OMHI] OMHI coefficients: ${OMHI_COEF}"
echo "[OMHI] OMHI universe: ${OMHI_UNIVERSE}"

# -----------------------------
# 4. Pre-flight checks
# -----------------------------
check_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "[ERROR] Command not found: $1"
    exit 1
  fi
}

check_file() {
  if [ ! -s "$1" ]; then
    echo "[ERROR] File not found or empty: $1"
    exit 1
  fi
}

check_command fastp
check_command bowtie2
check_command metaphlan
check_command merge_metaphlan_tables.py
check_command Rscript

check_file "${SAMPLE_SHEET}"
check_file "${OMHI_COEF}"
check_file "${OMHI_UNIVERSE}"
check_file "${LEVEL_EXTRACTOR_SCRIPT}"
check_file "${CALC_OMHI_SCRIPT}"

if [ ! -d "${MPA_DB}" ]; then
  echo "[ERROR] MetaPhlAn4 database directory does not exist: ${MPA_DB}"
  exit 1
fi

if [ ! -s "${HOST_INDEX}.1.bt2" ] && [ ! -s "${HOST_INDEX}.1.bt2l" ]; then
  echo "[ERROR] Bowtie2 host index prefix seems invalid: ${HOST_INDEX}"
  echo "[ERROR] Expected files like ${HOST_INDEX}.1.bt2 or ${HOST_INDEX}.1.bt2l"
  exit 1
fi

# Check sample sheet header
HEADER=$(head -n 1 "${SAMPLE_SHEET}" | tr -d '\r')
if [ "${HEADER}" != $'sample_id\tr1\tr2' ]; then
  echo "[ERROR] samples.tsv header must be exactly:"
  echo -e "sample_id\tr1\tr2"
  echo "[ERROR] Current header:"
  echo "${HEADER}"
  exit 1
fi

# -----------------------------
# 5. Process each sample
# -----------------------------
tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r SAMPLE_ID R1 R2
do
  SAMPLE_ID=$(echo "${SAMPLE_ID}" | tr -d '\r')
  R1=$(echo "${R1}" | tr -d '\r')
  R2=$(echo "${R2}" | tr -d '\r')

  if [ -z "${SAMPLE_ID}" ]; then
    continue
  fi

  echo ""
  echo "========================================================="
  echo "[OMHI] Processing sample: ${SAMPLE_ID}"
  echo "========================================================="

  check_file "${R1}"
  check_file "${R2}"

  CLEAN_R1="${QC_DIR}/${SAMPLE_ID}.clean.R1.fq.gz"
  CLEAN_R2="${QC_DIR}/${SAMPLE_ID}.clean.R2.fq.gz"

  NOHOST_R1="${NOHOST_DIR}/${SAMPLE_ID}.nohost.R1.fq.gz"
  NOHOST_R2="${NOHOST_DIR}/${SAMPLE_ID}.nohost.R2.fq.gz"

  MPA_PROFILE="${MPA_DIR}/${SAMPLE_ID}_profile.tsv"
  MPA_BOWTIE2_OUT="${MPA_DIR}/${SAMPLE_ID}.metaphlan.bowtie2.bz2"

  # -----------------------------
  # 5.1 QC with fastp
  # -----------------------------
  if [ -s "${CLEAN_R1}" ] && [ -s "${CLEAN_R2}" ]; then
    echo "[OMHI] QC output exists. Skipping fastp: ${SAMPLE_ID}"
  else
    echo "[OMHI] Running fastp: ${SAMPLE_ID}"

    fastp \
      -i "${R1}" \
      -I "${R2}" \
      -o "${CLEAN_R1}" \
      -O "${CLEAN_R2}" \
      --detect_adapter_for_pe \
      --thread "${THREADS}" \
      --html "${QC_DIR}/${SAMPLE_ID}.fastp.html" \
      --json "${QC_DIR}/${SAMPLE_ID}.fastp.json" \
      > "${LOG_DIR}/${SAMPLE_ID}.fastp.log" 2>&1
  fi

  # -----------------------------
  # 5.2 Remove host reads with Bowtie2
  # -----------------------------
  if [ -s "${NOHOST_R1}" ] && [ -s "${NOHOST_R2}" ]; then
    echo "[OMHI] No-host reads exist. Skipping Bowtie2 host removal: ${SAMPLE_ID}"
  else
    echo "[OMHI] Removing host reads with Bowtie2: ${SAMPLE_ID}"

    bowtie2 \
      -x "${HOST_INDEX}" \
      -1 "${CLEAN_R1}" \
      -2 "${CLEAN_R2}" \
      --very-sensitive \
      --threads "${THREADS}" \
      --un-conc-gz "${NOHOST_DIR}/${SAMPLE_ID}.nohost.R%.fq.gz" \
      -S /dev/null \
      > "${LOG_DIR}/${SAMPLE_ID}.bowtie2_host.log" 2>&1

    if [ ! -s "${NOHOST_R1}" ] || [ ! -s "${NOHOST_R2}" ]; then
      echo "[ERROR] Host-removed FASTQ files were not generated for sample: ${SAMPLE_ID}"
      echo "[ERROR] Expected:"
      echo "  ${NOHOST_R1}"
      echo "  ${NOHOST_R2}"
      exit 1
    fi
  fi

  # -----------------------------
  # 5.3 MetaPhlAn4 profiling
  # -----------------------------
  if [ -s "${MPA_PROFILE}" ]; then
    echo "[OMHI] MetaPhlAn4 profile exists. Skipping MetaPhlAn4: ${SAMPLE_ID}"
  else
    echo "[OMHI] Running MetaPhlAn4: ${SAMPLE_ID}"

    metaphlan \
      "${NOHOST_R1},${NOHOST_R2}" \
      --input_type fastq \
      --bowtie2db "${MPA_DB}" \
      --index "${MPA_INDEX}" \
      --nproc "${THREADS}" \
      --bowtie2out "${MPA_BOWTIE2_OUT}" \
      -o "${MPA_PROFILE}" \
      > "${LOG_DIR}/${SAMPLE_ID}.metaphlan4.log" 2>&1
  fi

done

# -----------------------------
# 6. Merge MetaPhlAn4 profiles
# -----------------------------
MERGED_PROFILE="${MERGED_DIR}/metaphlan4_merged.tsv"
SPECIES_TABLE="${MERGED_DIR}/taxonomy_species.tsv"
OMHI_OUTPUT="${OMHI_DIR}/omhi_scores.tsv"

echo ""
echo "========================================================="
echo "[OMHI] Merging MetaPhlAn4 profiles"
echo "========================================================="

PROFILE_COUNT=$(ls "${MPA_DIR}"/*_profile.tsv 2>/dev/null | wc -l)

if [ "${PROFILE_COUNT}" -eq 0 ]; then
  echo "[ERROR] No MetaPhlAn4 profile files found in: ${MPA_DIR}"
  exit 1
fi

merge_metaphlan_tables.py "${MPA_DIR}"/*_profile.tsv \
  > "${MERGED_PROFILE}"

check_file "${MERGED_PROFILE}"

# -----------------------------
# 7. Extract taxonomic levels from MetaPhlAn4 merged table
# -----------------------------
echo ""
echo "========================================================="
echo "[OMHI] Extracting taxonomic levels from MetaPhlAn4 table"
echo "========================================================="

LEVEL_DIR="${MERGED_DIR}/level"
mkdir -p "${LEVEL_DIR}"

python "${LEVEL_EXTRACTOR_SCRIPT}" \
  "${MERGED_PROFILE}" \
  "${LEVEL_DIR}"

echo "[OMHI] Files generated by level extractor:"
find "${LEVEL_DIR}" -maxdepth 1 -type f -print

# Try common species-level output names
if [ -s "${LEVEL_DIR}/species.tsv" ]; then
  SPECIES_TABLE="${LEVEL_DIR}/species.tsv"
elif [ -s "${LEVEL_DIR}/taxonomy_species.tsv" ]; then
  SPECIES_TABLE="${LEVEL_DIR}/taxonomy_species.tsv"
elif [ -s "${LEVEL_DIR}/metaphlan_species.tsv" ]; then
  SPECIES_TABLE="${LEVEL_DIR}/metaphlan_species.tsv"
else
  echo "[ERROR] Cannot find species-level table in: ${LEVEL_DIR}"
  echo "[ERROR] Expected one of:"
  echo "  ${LEVEL_DIR}/species.tsv"
  echo "  ${LEVEL_DIR}/taxonomy_species.tsv"
  echo "  ${LEVEL_DIR}/metaphlan_species.tsv"
  echo "[ERROR] Please check the output filename of metaphlan_level_extractor.py"
  exit 1
fi

check_file "${SPECIES_TABLE}"

echo "[OMHI] Species table used for OMHI:"
echo "  ${SPECIES_TABLE}"

# -----------------------------
# 8. Calculate OMHI
# -----------------------------
echo ""
echo "========================================================="
echo "[OMHI] Calculating OMHI scores"
echo "========================================================="

Rscript "${CALC_OMHI_SCRIPT}" \
  --input "${SPECIES_TABLE}" \
  --coef "${OMHI_COEF}" \
  --universe "${OMHI_UNIVERSE}" \
  --output "${OMHI_OUTPUT}"

check_file "${OMHI_OUTPUT}"

echo ""
echo "[OMHI] Pipeline finished at: $(date)"
echo "[OMHI] Final output:"
echo "  ${OMHI_OUTPUT}"
echo ""
echo "[OMHI] Preview:"
head "${OMHI_OUTPUT}"
