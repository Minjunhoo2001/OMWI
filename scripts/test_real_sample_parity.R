#!/usr/bin/env Rscript
# Run from repository root, without refitting.
suppressPackageStartupMessages(library(data.table))
for (split in c("train", "external")) {
  expected <- if (split == "train") "species_train_in_sample_scores.csv" else "external_independent_predictions.csv"
  out <- tempfile(fileext=".tsv")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("scripts/calc_omwi.R", "--input", shQuote(paste0("manuscript/data/abundance/", split, "_species.tsv")), "--output", shQuote(out)))
  if (status != 0L) stop("Scoring failed: ", split)
  a <- fread(out); b <- fread(file.path("manuscript/data/model_outputs/final", expected))
  stopifnot(!anyDuplicated(a$sample_id), !anyDuplicated(b$sample_id), setequal(a$sample_id,b$sample_id))
  i <- match(b$sample_id, a$sample_id)
  err <- max(abs(a$OMWI[i]-b$OMHI))
  if (!is.finite(err) || err > 1e-10) stop("Score mismatch: ",split,"; error=",err)
  cat(split, ": ", nrow(b), " matched; max absolute error=", format(err, digits=12), "\n", sep="")
  unlink(out)
}
