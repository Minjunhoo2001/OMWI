# Frozen OMWI model

The final 80-species model from the supplied compact Raw A output. Both coefficient TSVs include the intercept; they are identical because all 80 coefficients are nonzero. `omwi_species_universe.tsv` contains only the 80 panel taxa, in model order. `omwi_final_model.rds` is the original frozen RDS (a list with `fit`, `lambda`, `panel_taxa`, and other metadata); its internal training paths are historical, not runtime requirements for coefficient-based scoring.

For each sample, align to the 80 panel taxa, fill absent taxa with zero, divide by the sum across the 80 taxa (if positive), multiply by 100, add 1e-6, take natural logs, and subtract their mean. Score using the intercept and exported coefficients without additional standardization. An all-zero panel produces the intercept, but is flagged as unreliable. Do not interpret it as clinical disease. No MMUPHin is required.
