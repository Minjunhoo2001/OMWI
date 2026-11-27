#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for f in scripts/calc_omwi.R models/model_config.yml models/omwi_species_coefficients.tsv models/omwi_species_universe.tsv; do
  [[ -s "$f" ]] || { echo "Missing $f" >&2; exit 2; }
done
Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE))'
Rscript scripts/test_scoring.R
