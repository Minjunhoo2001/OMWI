#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

required_commands=(fastp bowtie2 metaphlan merge_metaphlan_tables.py python Rscript)
required_files=(
  "scripts/calc_omwi.R"
  "scripts/metaphlan_level_extractor.py"
  "models/omwi_species_coefficients.tsv"
  "models/omwi_species_coefficients_full.tsv"
  "models/omwi_species_model_config.tsv"
  "models/omwi_species_universe.tsv"
)

failed=0

for command_name in "${required_commands[@]}"; do
  if command -v "${command_name}" >/dev/null 2>&1; then
    printf '[OK] command: %s\n' "${command_name}"
  else
    printf '[MISSING] command: %s\n' "${command_name}" >&2
    failed=1
  fi
done

for file_name in "${required_files[@]}"; do
  if [[ -s "${file_name}" ]]; then
    printf '[OK] file: %s\n' "${file_name}"
  else
    printf '[MISSING] file: %s\n' "${file_name}" >&2
    failed=1
  fi
done

if [[ "${failed}" -ne 0 ]]; then
  printf '\nOMWI installation check failed. Resolve the missing items above.\n' >&2
  exit 1
fi

printf '\nOMWI installation check passed.\n'
printf 'External Bowtie2 and MetaPhlAn4 databases are validated when the full pipeline starts.\n'
