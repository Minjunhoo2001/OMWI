#!/usr/bin/env Rscript
# Run from repository root.
suppressPackageStartupMessages(library(data.table))
p <- tempfile(fileext=".tsv")
status <- system2(file.path(R.home("bin"),"Rscript"), c("scripts/calc_omwi.R","--input","example/example_species_abundance.tsv","--output",shQuote(p)))
if (status!=0L) stop("Scoring command failed")
a <- fread(p); b <- fread("example/expected_omwi_scores.tsv")
stopifnot(identical(a$sample_id,b$sample_id), max(abs(a$OMWI-b$OMWI))<1e-10, max(abs(a$prob_healthy-b$prob_healthy))<1e-10)
unlink(p)
cat("Synthetic score parity passed.\n")
