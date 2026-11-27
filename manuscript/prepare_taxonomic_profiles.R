#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

# 按每一行 Taxon 的末端层级提取，不将下级行重复汇总。
# 输入：Mixed 的两个已锁定数据集。无需重读原始 5954 样本 taxonomy。
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
SPLIT_DIR <- file.path(ROOT_DIR, "results/preparation/mixed_split")
OUT_DIR <- file.path(ROOT_DIR, "results/preparation/abundance")
LEVELS <- c("phylum", "class", "order", "family", "genus", "species", "SGB")
# 本版保留全部分类单元：不删除任何属/物种/SGB。
DROP_GENERA <- character(0)

read_otu <- function(file) {
  if (!file.exists(file)) stop("找不到文件：", file)
  x <- fread(file, check.names = FALSE)
  if (!"Taxon" %in% names(x)) stop("缺少 Taxon 列：", file)
  if (anyNA(x$Taxon) || any(trimws(x$Taxon) == "") || anyDuplicated(x$Taxon))
    stop("Taxon 存在缺失、空值或重复：", file)
  if (anyDuplicated(names(x))) stop("列名重复：", file)
  x[, Taxon := as.character(Taxon)]
  x
}
DROP_GENERA <- c(
  "Bradyrhizobium",
  "Ralstonia",
  "Sphingomonas",
  "Methylobacterium",
  "Rhizobium",
  "Mesorhizobium",
  "Ottowia"
)

remove_genera <- function(dt) {
  pat <- paste0(
    "(^|\\|)g__(",
    paste(DROP_GENERA, collapse = "|"),
    ")($|\\|)"
  )

  flag <- grepl(pat, dt$Taxon, perl = TRUE)

  cat("删除分类行数：", sum(flag), "\n")
  dt[!flag]
}
prefix_map <- c(phylum="p__", class="c__", order="o__", family="f__",
                genus="g__", species="s__", SGB="t__")

extract_terminal_level <- function(dt, level) {
  prefix <- prefix_map[[level]]
  pat <- paste0("(^|\\|)", prefix, "[^|]+$")
  idx <- grepl(pat, dt$Taxon, perl = TRUE)
  z <- copy(dt[idx])
  if (nrow(z) == 0L) stop("未找到末端层级为 ", level, " 的分类行，检查 Taxon 格式")
  # 末端标签之前的谱系全部去掉，再删前缀。
  z[, Taxon := sub(paste0("^.*(?:^|\\|)", prefix), "", Taxon, perl = TRUE)]
  if (anyNA(z$Taxon) || any(z$Taxon == "") || anyDuplicated(z$Taxon))
    stop(level, " 提取的名称缺失/重复：请检查同名分类单元，不自动相加")
  cols <- setdiff(names(z), "Taxon")
  for (cc in cols) {
    v <- suppressWarnings(as.numeric(z[[cc]]))
    if (anyNA(v) || any(!is.finite(v)) || any(v < 0))
      stop(level, ": 样本 ", cc, " 含 NA / 非数值 / 负丰度")
    set(z, j = cc, value = v)
  }
  setcolorder(z, c("Taxon", cols))
  z[]
}

# 不执行属过滤；保留所有原始分类行。
train <- read_otu(file.path(SPLIT_DIR, "otu_train_species.tsv"))
external <- read_otu(file.path(SPLIT_DIR, "otu_external_species.tsv"))
train <- remove_genera(train)
external <- remove_genera(external)
if (!identical(train$Taxon, external$Taxon))
  stop("Train 和 External 的 Taxon 列顺序不同，请检查拆分源文件")
if (length(intersect(names(train)[-1L], names(external)[-1L])) > 0L)
  stop("Train 和 External 样本 ID 有重叠")

stopifnot(identical(train$Taxon, external$Taxon))
all_dt <- cbind(copy(train), external[, -1L, with = FALSE])

# 修复上一版本 get('external') 找不到对象的问题：显式使用数据列表。
datasets <- list(train = train, external = external, all = all_dt)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
summary <- list()
for (lv in LEVELS) {
  for (nm in names(datasets)) {
    z <- extract_terminal_level(datasets[[nm]], lv)
    f <- file.path(OUT_DIR, paste0(nm, "_", lv, ".tsv"))
    fwrite(z, f, sep = "\t", quote = FALSE)
    summary[[length(summary) + 1L]] <- data.table(
      dataset = nm, level = lv, n_taxa = nrow(z), n_samples = ncol(z) - 1L
    )
    cat(sprintf("%s: %d taxa × %d samples -> %s\n",
                paste(nm, lv, sep = "_"), nrow(z), ncol(z) - 1L, f))
  }
}
summary_dt <- rbindlist(summary)
expected <- c(phylum=25L, class=288L, order=323L, family=389L, genus=1183L, species=2523L, SGB=2711L)
cat("\n注意：预计类别数以当前审计输出为参考，若 taxonomy 更新可不一致。\n")
print(summary_dt[dataset=="train", .(level,n_taxa,n_samples)])
fwrite(summary_dt, file.path(OUT_DIR, "level_extraction_summary.csv"))
if (nrow(summary_dt) != length(LEVELS) * length(datasets))
  stop("各层级输出文件数量不完整")
cat("\n完成！7 个层级 × 3 个数据集 = 21 个 TSV。目录：", OUT_DIR, "\n")
