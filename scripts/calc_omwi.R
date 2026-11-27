#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i+1L], "--")) stop("Missing value for ", flag)
  args[i+1L]
}
if ("--help" %in% args || is.null(get_arg("--input"))) {
  cat("Usage: Rscript scripts/calc_omwi.R --input species.tsv [--coef models/omwi_species_coefficients.tsv] [--universe models/omwi_species_universe.tsv] [--output omwi_scores.tsv]\n")
  quit(status=if ("--help" %in% args) 0L else 2L)
}
input <- get_arg("--input")
coef_file <- get_arg("--coef", "models/omwi_species_coefficients.tsv")
universe_file <- get_arg("--universe", "models/omwi_species_universe.tsv")
output <- get_arg("--output", "omwi_scores.tsv")
pc <- as.numeric(get_arg("--pseudocount", "1e-6"))
if (!is.finite(pc) || pc != 1e-6) stop("Frozen model requires pseudocount 1e-6 on percentage scale")
d <- fread(input, data.table=FALSE, check.names=FALSE)
if (ncol(d)<2L || nrow(d)<1L) stop("Expected species rows and sample columns")
ids <- trimws(names(d)[-1L])
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid or duplicate sample IDs")
f <- trimws(as.character(d[[1L]]))
if (anyNA(f) || any(!nzchar(f))) stop("Invalid feature names")
# Accept bare species IDs or terminal species-level MetaPhlAn clades.
# Never sum strain/SGB descendants or higher ranks into species-level rows.
keep <- grepl("(^|\\|)s__[^|]+$", f) | (!grepl("\\|", f) & !grepl("^[dkpcofgst]__", f))
f <- sub("^.*\\|s__", "", f[keep]); f <- sub("^s__", "", f)
if (!length(f)) stop("No species-level rows found")
a <- vapply(d[keep,-1L,drop=FALSE], function(v) suppressWarnings(as.numeric(v)), numeric(sum(keep)))
a <- matrix(a,nrow=sum(keep),ncol=length(ids),dimnames=list(f,ids))
if (any(!is.finite(a)) || any(a<0)) stop("Abundance values must be finite and nonnegative")
if (anyDuplicated(f)) a <- rowsum(a,group=f,reorder=FALSE)
u <- fread(universe_file,data.table=FALSE)
if (!"feature" %in% names(u)) stop("Universe file must contain feature column")
panel <- as.character(u$feature)
if (length(panel)!=80L || anyNA(panel) || any(!nzchar(panel)) || anyDuplicated(panel)) stop("Expected fixed 80-species universe")
cf <- fread(coef_file,data.table=FALSE)
if (!all(c("feature","coef") %in% names(cf))) stop("Coefficient file needs feature and coef columns")
cf$feature <- as.character(cf$feature); cf$coef <- as.numeric(cf$coef)
if (nrow(cf)!=81L || anyNA(cf$feature) || anyDuplicated(cf$feature) || any(!is.finite(cf$coef)) || !setequal(cf$feature,c("(Intercept)",panel))) stop("Coefficient file and universe do not match")
beta <- cf$coef[match(panel,cf$feature)]
if (sum(beta>0)!=39L || sum(beta<0)!=41L) stop("Expected 39 positive and 41 negative species coefficients")
intercept <- cf$coef[match("(Intercept)",cf$feature)]
matched <- intersect(panel,rownames(a))
if (!length(matched)) stop("No panel taxa matched; check taxonomy identifiers/database")
x <- matrix(0,nrow=80L,ncol=length(ids),dimnames=list(panel,ids))
x[matched,] <- a[matched,,drop=FALSE]
s <- colSums(x)
if (any(!is.finite(s))) stop("Non-finite panel sums")
good <- s>0
if (any(good)) x[,good] <- sweep(x[,good,drop=FALSE],2L,s[good],"/")
lx <- log(x*100+pc)
z <- sweep(lx,2L,colMeans(lx),"-")
score <- intercept+as.numeric(crossprod(beta,z))
# Preserve the training implementation for all-zero panel profiles.
score[!good] <- intercept
out <- data.frame(sample_id=ids,OMWI=score,
  predicted_state=ifelse(score>0,"health_associated","non_healthy_associated"),
  matched_universe_taxa=length(matched),total_universe_taxa=length(panel),
  prob_healthy=plogis(score))
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
fwrite(out,output,sep="\t")
cat("[OMWI] Wrote ",nrow(out)," scores to ",output,"\n",sep="")
