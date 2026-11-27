#!/usr/bin/env bash
# Paired-end FASTQ convenience workflow for the frozen OMWI model.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ $# -eq 2 ]] || { echo 'Usage: bash workflow/run_omwi_pipeline.sh samples.tsv NEW_OUTPUT_DIR' >&2; exit 2; }
SHEET="$1"; OUT="$2"
: "${HOST_INDEX:?Set HOST_INDEX to the GRCh38 Bowtie2 index prefix}"
: "${MPA_DB:?Set MPA_DB to the MetaPhlAn database directory}"
INDEX=mpa_vJun23_CHOCOPhlAnSGB_202403
THREADS="${THREADS:-8}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'THREADS must be a positive integer' >&2; exit 2; }
for c in fastp bowtie2 samtools metaphlan merge_metaphlan_tables.py python Rscript; do command -v "$c" >/dev/null || { echo "Missing $c" >&2; exit 2; }; done
[[ -s "${HOST_INDEX}.1.bt2" || -s "${HOST_INDEX}.1.bt2l" ]] || { echo 'Invalid HOST_INDEX prefix' >&2; exit 2; }
[[ -s "$MPA_DB/$INDEX.pkl" ]] || { echo "Missing $INDEX.pkl in MPA_DB" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo 'Choose a new output directory to avoid mixing profiles from different runs' >&2; exit 2; }
# Validate the complete sample sheet before starting any processing.
python - "$SHEET" <<'PY'
import csv, os, re, sys
with open(sys.argv[1], newline='', encoding='utf-8-sig') as f:
    rows = list(csv.reader(f, delimiter='\t'))
if not rows or rows[0] != ['sample_id','r1','r2'] or len(rows)<2:
    raise SystemExit('Expected sample_id, r1, r2 tab-delimited header and at least one sample')
seen=set()
for row in rows[1:]:
    if len(row)!=3 or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',row[0]) or row[0] in seen:
        raise SystemExit(f'Invalid or duplicate sample row: {row!r}')
    if any(not os.path.isfile(p) or not os.path.getsize(p) or ',' in p for p in row[1:]):
        raise SystemExit(f'Missing/empty FASTQ or unsupported comma in path: {row!r}')
    if os.path.realpath(row[1])==os.path.realpath(row[2]):
        raise SystemExit('R1 and R2 must be distinct files')
    seen.add(row[0])
PY
bash "$ROOT/workflow/test_installation.sh"
mkdir -p "$OUT"/{01.qc,02.nohost,03.metaphlan,04.merged,05.omwi,logs}
profiles=()
while IFS=$'\t' read -r id r1 r2 || [[ -n "${id:-}" ]]; do
  r2="${r2%$'\r'}"
  qc="$OUT/01.qc/$id"; nh="$OUT/02.nohost/$id"; mp="$OUT/03.metaphlan/$id"
  fastp -i "$r1" -I "$r2" -o "$qc.R1.fq.gz" -O "$qc.R2.fq.gz" \
    --detect_adapter_for_pe --thread "$THREADS" --html "$qc.html" --json "$qc.json" > "$OUT/logs/$id.fastp.log" 2>&1
  # Retain only primary read pairs for which both mates are unmapped.
  bowtie2 --very-sensitive -x "$HOST_INDEX" -1 "$qc.R1.fq.gz" -2 "$qc.R2.fq.gz" -p "$THREADS" 2> "$OUT/logs/$id.host.log" | \
    samtools view -u -f 12 -F 2304 - | samtools sort -n -@ "$THREADS" -o "$nh.bam" -
  samtools fastq -@ "$THREADS" -n -1 "$nh.R1.fq.gz" -2 "$nh.R2.fq.gz" -0 /dev/null -s /dev/null "$nh.bam" 2>> "$OUT/logs/$id.host.log"
  rm "$nh.bam"
  metaphlan "$nh.R1.fq.gz,$nh.R2.fq.gz" --input_type fastq --bowtie2db "$MPA_DB" --index "$INDEX" \
    --nproc "$THREADS" --bowtie2out "$mp.bowtie2.bz2" -o "$mp.tsv" > "$OUT/logs/$id.metaphlan.log" 2>&1
  profiles+=("$mp.tsv")
done < <(tail -n +2 "$SHEET")
merged="$OUT/04.merged/metaphlan4_merged.tsv"
merge_metaphlan_tables.py "${profiles[@]}" > "$merged"
python "$ROOT/scripts/metaphlan_level_extractor.py" "$merged" "$OUT/04.merged"
Rscript "$ROOT/scripts/calc_omwi.R" --input "$OUT/04.merged/metaphlan4_merged_species.tsv" \
  --coef "$ROOT/models/omwi_species_coefficients.tsv" --universe "$ROOT/models/omwi_species_universe.tsv" \
  --output "$OUT/05.omwi/omwi_scores.tsv"
