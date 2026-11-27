#!/usr/bin/env Rscript

# =========================================================
# Supplementary Fig. 3 | OMWI differences and discrimination by age, sex, geography and oral sampling site
# Twelve panels (4 rows x 3 columns) using the ORIGINAL figure style:
# a age-group OMWI; b health status within age; c age-specific ROC
# d sex-group OMWI; e health status within sex; f sex-specific ROC
# g continent OMWI; h health status within continent; i continent ROC
# j site OMWI; k health status within site; l site ROC
# Panels g/j: BH-adjusted pairwise Wilcoxon significance stars via ggpubr; no overall P annotation.
# Site panels use Saliva, Tongue dorsum, Pharyngeal swab (oral rinse excluded only from site panels).
# Geographic panels use all available model-score samples FROM THREE CONTINENTS
# (Asia, Europe, Americas); Oceania is excluded ONLY from geographic panels.
# All source files are read from ROOT_DIR (no old Oral_proxy2.0 inputs).
# Source model: data/model_outputs / final
# NOTE: training scores are in-sample. Pooled ROC is descriptive,
# not an independent validation estimate. Set ANALYSIS_SCOPE="test"
# for independent validation-only AGE/SEX panels; geography and site stay full-sample.
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(pROC)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. Paths and analysis parameters
# =========================================================
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
OUT_DIR  <- file.path(ROOT_DIR, "Supplementaryfigure", "bias")

# is present under ROOT_DIR, the two input CSVs are extracted to OUT_DIR.
MODEL_DIR_NAME <- "data/model_outputs"
# Geographic and sampling-site metadata are used as lookup tables, NOT
# as the source of OMWI scores. The metadata may contain MORE samples than
# the score CSVs, whose true training row count must be checked at load time.
# Every scored sample must have a unique metadata match.
# The existing path below can be used if all model IDs exist in this file;
# its folder name alone does not prove that the contents are outdated.
GEO_META_FILE <- file.path(ROOT_DIR, "data/metadata", "mix_sited_metadata.csv")

# "all" = new-model training (in-sample) + independent external samples;
# "test" = external-only, recommended for independent validation claims.
ANALYSIS_SCOPE <- "all"
SHOW_INTERACTION_ON_ROC <- FALSE   # reference PDF: only ROC/AUC legend

MIN_PER_CLASS_FOR_VIOLIN <- 5L
MIN_PER_CLASS_FOR_AUC    <- 10L
ROC_BOOT_N               <- 2000L
SET_SEED                  <- 704L

if (!dir.exists(ROOT_DIR)) {
  stop("ROOT_DIR does not exist: ", ROOT_DIR)
}
for (subdir in c("", "individual_panels", "tables")) {
  dir.create(file.path(OUT_DIR, subdir), recursive = TRUE, showWarnings = FALSE)
}

# =========================================================
# 1. Colors and plotting theme, matched to Fig. 4 (larger publication labels)
# =========================================================
COL_1 <- "#4DBBD5"
COL_2 <- "#FEE5D9"
COL_3 <- "#FCAE91"
COL_4 <- "#FB6A4A"
COL_5 <- "#CB181D"

COL_HEALTH  <- COL_1
COL_DISEASE <- COL_5
# Blue/red are reserved exclusively for Healthy vs Disease in panels b/e/h.
# All age/sex/continent grouping colors (a/d/g and their matching ROC lines)
# use sage green, ochre, and lilac/purple to distinguish the comparisons.
AGE_COLORS <- c(
  "18–44" = "#669B7A",  # sage green
  "45–59" = "#D7A650",  # ochre
  "≥60"   = "#9879B5"   # muted purple
)

SEX_COLORS <- c(
  "Male"   = "#6CA88E",  # muted green
  "Female" = "#B397CF"   # soft lilac
)

HD_COLORS <- c(
  "Healthy" = COL_HEALTH,
  "Disease" = COL_DISEASE
)

GEO_LEVELS <- c("Asia", "Europe", "Americas")
GEO_COLORS <- c(
  "Asia"     = "#9076B2",  # purple
  "Europe"   = "#69A487",  # muted green
  "Americas" = "#D0A04F"   # ochre
)

# Text sizes are in points except ggplot geom_text size (millimetres).
AXIS_FONT_SIZE      <- 18
STRIP_FONT_SIZE     <- 18
LEGEND_FONT_SIZE    <- 10
ANNOTATION_SIZE_MM  <- 5.0
PANEL_TAG_FONT_SIZE <- 18

# Country of sampling, NOT ethnicity. Unknown labels are audited and rejected.
GEO_COUNTRY_MAP <- c(
  "China" = "Asia", "Japan" = "Asia", "Singapore" = "Asia",
  "Philippines" = "Asia", "South Korea" = "Asia", "India" = "Asia",
  "Germany" = "Europe", "Denmark" = "Europe", "Spain" = "Europe",
  "Italy" = "Europe", "France" = "Europe", "Luxembourg" = "Europe",
  "Netherlands" = "Europe", "United Kingdom" = "Europe",
  "Sweden" = "Europe", "Austria" = "Europe", "Finland" = "Europe",
  "United States" = "Americas", "USA" = "Americas",
  "United States of America" = "Americas", "Canada" = "Americas",
  "Peru" = "Americas", "Brazil" = "Americas",
  "Fiji" = "Oceania", "Australia" = "Oceania", "New Zealand" = "Oceania"
)

theme_cell <- function(base_size = AXIS_FONT_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = AXIS_FONT_SIZE),
      axis.title = element_text(color = "black", size = AXIS_FONT_SIZE),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", color = "black", size = STRIP_FONT_SIZE),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", color = "black", size = LEGEND_FONT_SIZE),
      legend.text = element_text(color = "black", size = LEGEND_FONT_SIZE),
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5),
      plot.tag = element_text(face = "bold", color = "black", size = PANEL_TAG_FONT_SIZE),
      plot.tag.position = c(0.01, 0.99)
    )
}

format_p_label <- function(p, prefix = "P") {
  if (length(p) == 0L || !is.finite(p[1L])) return(paste0(prefix, " = NA"))
  p <- p[1L]
  if (p < 0.001) return(paste0(prefix, " = ", sprintf("%.2E", p)))
  paste0(prefix, " = ", sprintf("%.3f", p))
}

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

# Convert mixed UTF-8 / GBK / GB18030 metadata text to valid UTF-8.
# Conversion is done element by element because one metadata column may
# contain strings written with different encodings.
safe_text <- function(x) {
  x_chr <- as.character(x)
  out <- rep(NA_character_, length(x_chr))
  
  for (i in seq_along(x_chr)) {
    xi <- x_chr[i]
    if (is.na(xi)) next
    
    is_valid_utf8 <- tryCatch(
      isTRUE(validUTF8(xi)),
      error = function(e) FALSE
    )
    
    if (is_valid_utf8) {
      yi <- enc2utf8(xi)
    } else {
      yi <- suppressWarnings(
        iconv(xi, from = "GB18030", to = "UTF-8", sub = "byte")
      )
      
      if (is.na(yi)) {
        yi <- suppressWarnings(
          iconv(xi, from = "CP936", to = "UTF-8", sub = "byte")
        )
      }
      
      if (is.na(yi)) {
        yi <- suppressWarnings(
          iconv(xi, from = "latin1", to = "UTF-8", sub = "byte")
        )
      }
    }
    
    out[i] <- yi
  }
  
  out <- gsub("\u00A0", " ", out, fixed = TRUE)
  trimws(out)
}

# Strict numeric conversion after repairing text encoding.
# Values such as ">=18", "≥18", "18-44" or "unknown" are intentionally
# returned as NA because they are not exact individual-level ages.
safe_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  
  x_fix <- safe_text(x)
  x_fix <- gsub(",", "", x_fix, fixed = TRUE)
  x_fix <- gsub("[[:space:]]+", "", x_fix, perl = TRUE)
  x_fix <- gsub("(years?|yrs?|yearold|岁)$", "", x_fix,
                ignore.case = TRUE, perl = TRUE)
  
  missing_text <- is.na(x_fix) | x_fix == "" |
    tolower(x_fix) %in% c("na", "n/a", "nan", "null", "none", "unknown", "missing")
  
  numeric_pattern <- "^[+-]?(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))([eE][+-]?[0-9]+)?$"
  is_exact_numeric <- !missing_text & grepl(numeric_pattern, x_fix, perl = TRUE)
  
  out <- rep(NA_real_, length(x_fix))
  out[is_exact_numeric] <- suppressWarnings(as.numeric(x_fix[is_exact_numeric]))
  out
}

pick_first_existing_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# =========================================================
# 2. Read scores and metadata from the new Oral v3.0 model only
# =========================================================
# The ZIP's final CSVs ALREADY contain sample_id, age_num, gender,
# disease, study_name and OMWI. Joining old 01.split64/04.lasso64 data
# would silently mix models and is therefore not performed.
resolve_input_csvs <- function() {
  requested <- c(
    "species_train_in_sample_scores.csv",
    "external_independent_predictions.csv"
  )
  
  # Option 1: directly extracted folder, e.g. ROOT_DIR/05.raw_A.../final/
  candidate_folders <- list.dirs(ROOT_DIR, recursive = FALSE, full.names = TRUE)
  candidate_folders <- candidate_folders[
    grepl("^05[.]raw_A_stability_subset_v1", basename(candidate_folders), ignore.case = TRUE)
  ]
  candidates <- unique(c(
    file.path(ROOT_DIR, MODEL_DIR_NAME, "final"),
    file.path(candidate_folders, "final"),
    file.path(ROOT_DIR, "final")
  ))
  for (d in candidates) {
    fs <- file.path(d, requested)
    if (all(file.exists(fs))) {
      message("Reading extracted model outputs in: ", d)
      return(stats::setNames(fs, requested))
    }
  }
  
  # Option 2: ZIP archive saved under ROOT_DIR
  zip_files <- list.files(ROOT_DIR, pattern = "^05[.]raw_A_stability_subset_v1.*[.]zip$",
                          full.names = TRUE, ignore.case = TRUE)
  if (length(zip_files) == 0L) {
    stop("Missing v3 score files. Extract data/model_outputs.zip ",
         "into ROOT_DIR, or put that ZIP directly in ROOT_DIR. ",
         "Required files: ", paste(requested, collapse = ", "))
  }
  zip_file <- zip_files[1L]
  archive_names <- utils::unzip(zip_file, list = TRUE)$Name
  members <- vapply(requested, function(nm) {
    hits <- archive_names[grepl(paste0("(^|/)final/", nm, "$"), archive_names)]
    if (length(hits) != 1L) {
      stop("Expected exactly one ", nm, " under final/ inside ", zip_file)
    }
    hits[1L]
  }, character(1))
  
  cache_dir <- file.path(OUT_DIR, "input_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_file, files = unname(members), exdir = cache_dir,
               overwrite = TRUE)
  files <- file.path(cache_dir, members)
  if (!all(file.exists(files))) {
    stop("ZIP extraction failed. Expected: ", paste(files, collapse = "; "))
  }
  message("Reading v3 model scores extracted from: ", zip_file)
  stats::setNames(files, requested)
}

# fread reads CSV bytes without assuming all rows use GB18030. Forcing a
# single decoding on mixed-encoding study-name fields can create an
# "invalid input" warning and even truncate the read. Text used in analyses
# is normalized with safe_text() below; numeric scores are never changed.
read_v3_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  ans <- tryCatch(
    data.table::fread(path, encoding = "unknown", data.table = FALSE,
                      check.names = FALSE, showProgress = FALSE),
    error = function(e) {
      stop("Could not parse score CSV: ", path, "\n", conditionMessage(e))
    }
  )
  if (nrow(ans) == 0L) stop("Score CSV contains no data rows: ", path)
  ans
}

input_files <- resolve_input_csvs()
train_df <- read_v3_csv(input_files[["species_train_in_sample_scores.csv"]])
test_df  <- read_v3_csv(input_files[["external_independent_predictions.csv"]])

required_cols <- c("sample_id", "OMHI", "disease", "study_name", "gender", "age_num")
# Repair mixed-encoding character columns and normalize IDs BEFORE checking
# uniqueness or overlap. This does not modify model predictions.
for (nm in c("sample_id", "disease", "study_name", "gender")) {
  if (nm %in% names(train_df)) train_df[[nm]] <- safe_text(train_df[[nm]])
  if (nm %in% names(test_df))  test_df[[nm]]  <- safe_text(test_df[[nm]])
}
if ("sample_id" %in% names(train_df)) train_df$sample_id <- trimws(train_df$sample_id)
if ("sample_id" %in% names(test_df))  test_df$sample_id  <- trimws(test_df$sample_id)
for (label in c("train", "test")) {
  x <- if (label == "train") train_df else test_df
  miss <- setdiff(required_cols, names(x))
  if (length(miss)) stop(label, " file missing columns: ", paste(miss, collapse = ", "))
  if (anyDuplicated(x$sample_id)) stop(label, " file contains duplicate sample_id values.")
}
if (length(intersect(as.character(train_df$sample_id), as.character(test_df$sample_id))) > 0) {
  stop("Train and external test files overlap in sample_id. Refusing to combine.")
}
# Keep only the analysis columns: other v3 metadata fields have different
# inferred types between train and test (e.g. Disease activity).
train_df <- train_df[, required_cols, drop = FALSE]
test_df  <- test_df[, required_cols, drop = FALSE]
train_df$dataset_split <- "train"
test_df$dataset_split  <- "test"
dat <- dplyr::bind_rows(train_df, test_df) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    OMWI = safe_numeric(OMHI),
    dataset_split = as.character(dataset_split)
  )
if (all(is.na(dat$OMWI))) stop("OMWI is entirely missing in v3 model outputs.")
message("[v3 input] train n=", nrow(train_df), "; external n=", nrow(test_df),
        "; total n=", nrow(dat))
message("[v3 input] The training-score CSV supplies ", nrow(train_df),
        " scored training samples. Any other training subjects present only in",
        " the metadata are NOT silently assigned OMWI scores.")
if (nrow(train_df) == 330L) {
  warning("[CHECK TRAINING CSV] Only 330 training rows were parsed. The earlier ",
          "read.csv(GB18030) produced an invalid-input warning and may have ",
          "stopped early. fread no longer forces that encoding, but verify ",
          "whether 330 is the intended FULL training-score count.",
          call. = FALSE)
}

# =========================================================
# 3. Standardize outcome, age, sex and study variables
# =========================================================
if (!"disease" %in% colnames(dat)) {
  stop("The metadata do not contain the disease column.")
}

raw_disease <- tolower(safe_text(dat$disease))
dat$health_status <- dplyr::case_when(
  raw_disease %in% c("healthy", "health", "control", "ctrl", "0") ~ "Healthy",
  raw_disease %in% c("disease", "diseased", "nonhealthy", "non-healthy", "case", "1") ~ "Disease",
  TRUE ~ NA_character_
)
dat$health_status <- factor(dat$health_status, levels = c("Healthy", "Disease"))
dat$disease01 <- ifelse(dat$health_status == "Disease", 1, ifelse(dat$health_status == "Healthy", 0, NA))

age_col <- pick_first_existing_col(dat, c("age_num", "age", "Age", "AGE", "age_years", "Age_years"))
if (is.na(age_col)) stop("No age column was found.")

age_text_fixed <- safe_text(dat[[age_col]])
dat$age_num <- safe_numeric(age_text_fixed)

age_invalid <- unique(age_text_fixed[
  !is.na(age_text_fixed) & age_text_fixed != "" & is.na(dat$age_num)
])
if (length(age_invalid) > 0) {
  message(
    "Age entries excluded because they are not exact numeric ages (n = ",
    sum(!is.na(age_text_fixed) & age_text_fixed != "" & is.na(dat$age_num)),
    "). Examples: ",
    paste(utils::head(age_invalid, 10), collapse = "; ")
  )
}
dat$age_group <- dplyr::case_when(
  !is.na(dat$age_num) & dat$age_num >= 18 & dat$age_num <= 44 ~ "18–44",
  !is.na(dat$age_num) & dat$age_num >= 45 & dat$age_num <= 59 ~ "45–59",
  !is.na(dat$age_num) & dat$age_num >= 60 ~ "≥60",
  TRUE ~ NA_character_
)
dat$age_group <- factor(dat$age_group, levels = c("18–44", "45–59", "≥60"), ordered = TRUE)
dat$age_order <- as.numeric(dat$age_group)

sex_col <- pick_first_existing_col(dat, c("gender", "Gender", "sex", "Sex"))
if (is.na(sex_col)) stop("No sex or gender column was found.")

raw_sex <- tolower(safe_text(dat[[sex_col]]))
dat$sex_group <- dplyr::case_when(
  raw_sex %in% c("male", "m", "man", "boy", "1", "男") ~ "Male",
  raw_sex %in% c("female", "f", "woman", "girl", "2", "女") ~ "Female",
  TRUE ~ NA_character_
)
dat$sex_group <- factor(dat$sex_group, levels = c("Male", "Female"))

study_col <- pick_first_existing_col(
  dat,
  c("study_name", "Study", "study", "study_id", "StudyID", "project", "Project", "cohort", "Cohort", "PMID")
)

if (!is.na(study_col)) {
  dat$study_use <- factor(as.character(dat[[study_col]]))
} else {
  dat$study_use <- factor("single_study")
}

if (ANALYSIS_SCOPE == "test") {
  dat_use <- dat %>% filter(dataset_split == "test")
} else if (ANALYSIS_SCOPE == "all") {
  dat_use <- dat
} else {
  stop("ANALYSIS_SCOPE must be either 'all' or 'test'.")
}

dat_use <- dat_use %>%
  filter(!is.na(OMWI), !is.na(health_status))

cat("\nAnalysis scope:", ANALYSIS_SCOPE, "\n")
cat("Samples with OMWI and outcome:", nrow(dat_use), "\n")

# =========================================================
# 3.1. Geographic metadata (full-sample analysis; does NOT affect age/sex)
# =========================================================
# A metadata file may describe the entire eligible cohort (e.g. 1793 rows),
# while the number of model scores must be established by the corrected CSV
# reader. The VALID matching rule is "model IDs are a subset of metadata IDs", not
# "all metadata IDs must have scores". Preserve the rest as an exclusion audit.
find_current_geo_metadata <- function(model_ids) {
  active_source_files <- vapply(sys.frames(), function(fr) {
    if (is.null(fr$ofile)) "" else as.character(fr$ofile)
  }, character(1))
  active_source_files <- active_source_files[nzchar(active_source_files)]
  source_dirs <- if (length(active_source_files)) {
    dirname(normalizePath(active_source_files, winslash = "/", mustWork = FALSE))
  } else character(0)
  
  user_home <- path.expand("~")
  win_home <- Sys.getenv("USERPROFILE", unset = "")
  candidate_dirs <- unique(c(
    ROOT_DIR,
    source_dirs,
    
    
    file.path(ROOT_DIR, "Supplementaryfigure", "bias"),
    file.path(user_home, "Downloads"),
    if (nzchar(win_home)) file.path(win_home, "Downloads") else character(0)
  ))
  candidate_names <- c("mix_sited_metadata.csv", "meta_oral_proxy_all(2).csv")
  candidates <- unique(c(GEO_META_FILE, unlist(lapply(candidate_dirs, function(d) {
    file.path(d, candidate_names)
  }), use.names = FALSE)))
  candidates <- candidates[file.exists(candidates)]
  
  load_and_validate <- function(path) {
    x <- tryCatch(
      data.table::fread(path, encoding = "unknown", data.table = FALSE,
                        check.names = FALSE, showProgress = FALSE),
      error = function(e) {
        message("[metadata] CSV parse failed: ", conditionMessage(e))
        NULL
      }
    )
    required <- c("sample_id", "country", "study_name", "disease", "body_site")
    if (is.null(x)) return(NULL)
    missing_columns <- setdiff(required, names(x))
    if (length(missing_columns)) {
      message("[metadata] Missing columns in ", path, ": ",
              paste(missing_columns, collapse = ", "))
      return(NULL)
    }
    if ("screening_status" %in% names(x)) {
      included <- !is.na(x$screening_status) &
        toupper(trimws(safe_text(x$screening_status))) == "INCLUDE"
      x <- x[included, , drop = FALSE]
    }
    x$sample_id <- trimws(safe_text(x$sample_id))
    if (anyNA(x$sample_id) || any(!nzchar(x$sample_id))) {
      message("[metadata] Missing/empty sample_id in ", path)
      return(NULL)
    }
    if (anyDuplicated(x$sample_id)) {
      message("[metadata] Duplicate sample_id in ", path)
      return(NULL)
    }
    missing_ids <- setdiff(model_ids, x$sample_id)
    if (length(missing_ids)) {
      message("[metadata] ", length(missing_ids), " scored sample IDs absent from ",
              path, "; first IDs: ",
              paste(head(missing_ids, 8L), collapse = ", "))
      return(NULL)
    }
    extra_ids <- setdiff(x$sample_id, model_ids)
    # Only rows with model scores enter downstream metadata joins. Matching
    # in model order avoids accidental row multiplication and label changes.
    matched <- x[match(model_ids, x$sample_id), , drop = FALSE]
    if (nrow(matched) != length(model_ids) ||
        !identical(as.character(matched$sample_id), as.character(model_ids))) {
      stop("Internal metadata ordering/ID match failed for ", path)
    }
    list(data = matched, n_metadata = nrow(x), n_extra = length(extra_ids),
         extra_data = x[x$sample_id %in% extra_ids,
                        c("sample_id", "study_name", "country", "body_site"),
                        drop = FALSE])
  }
  
  finalize_selection <- function(path, candidate) {
    message("[metadata] Using: ", normalizePath(path, winslash = "/"))
    message("[metadata] after screening: ", candidate$n_metadata,
            "; model IDs matched: ", nrow(candidate$data),
            "; metadata-only IDs (without scores): ", candidate$n_extra)
    data.table::fwrite(candidate$extra_data,
                       file.path(OUT_DIR, "tables", "metadata_only_unscored_samples_CHECK.csv"))
    data.table::fwrite(data.frame(sample_id = model_ids),
                       file.path(OUT_DIR, "tables", "model_ids_matched_to_metadata_CHECK.csv"))
    list(path = path, data = candidate$data,
         n_metadata = candidate$n_metadata, n_extra = candidate$n_extra)
  }
  
  for (path in candidates) {
    candidate <- load_and_validate(path)
    if (!is.null(candidate)) return(finalize_selection(path, candidate))
    message("[metadata] This candidate is not usable: ", path)
  }
  if (!interactive()) {
    stop("No metadata CSV contained all ", length(model_ids),
         " scored IDs with the required columns (sample_id, country,",
         " study_name, disease, body_site). Check candidate metadata files.")
  }
  message("[metadata] Select a CSV containing ALL ", length(model_ids),
          " scored IDs. Extra metadata-only samples are allowed.")
  chosen_path <- tryCatch(file.choose(), error = function(e) "")
  if (!nzchar(chosen_path)) stop("No metadata file selected.")
  chosen <- load_and_validate(chosen_path)
  if (is.null(chosen)) {
    stop("Selected metadata is missing scored IDs or required columns,",
         " or has duplicate IDs. Its TOTAL row count need not be ",
         length(model_ids), ".")
  }
  finalize_selection(chosen_path, chosen)
}

model_ids <- trimws(safe_text(dat$sample_id))
if (anyNA(model_ids) || any(!nzchar(model_ids)) || anyDuplicated(model_ids)) {
  stop("Model predictions have missing, empty, or duplicate sample IDs; do not join metadata.")
}
geo_selected <- find_current_geo_metadata(model_ids)
GEO_META_FILE <- geo_selected$path
geo_meta <- geo_selected$data
message("[metadata] matched rows entering geographic and site analyses: ",
        nrow(geo_meta), "; full metadata rows: ", geo_selected$n_metadata)
geo_meta <- geo_meta %>%
  mutate(
    sample_id = trimws(as.character(sample_id)),
    country = trimws(safe_text(country)),
    geo_study_name = trimws(safe_text(study_name)),
    geo_disease = tolower(trimws(safe_text(disease)))
  )

unknown_countries <- sort(unique(geo_meta$country[
  is.na(geo_meta$country) | is.na(GEO_COUNTRY_MAP[geo_meta$country])
]))
if (length(unknown_countries)) {
  fwrite(
    geo_meta %>% filter(is.na(country) | !country %in% names(GEO_COUNTRY_MAP)) %>%
      select(sample_id, geo_study_name, country),
    file.path(OUT_DIR, "tables", "geography_unmapped_countries_CHECK.csv")
  )
  stop("Unmapped country labels in current metadata: ",
       paste(unknown_countries, collapse = "; "),
       ". Verify sampling location before plotting.")
}
geo_meta$continent <- unname(GEO_COUNTRY_MAP[geo_meta$country])

# OMWI scores and outcomes always come from the frozen model, not metadata.
dat_geo <- dat %>%
  mutate(sample_id = trimws(as.character(sample_id))) %>%
  left_join(
    geo_meta %>% select(sample_id, country, geo_study_name, geo_disease, continent),
    by = "sample_id"
  )
missing_geo <- dat_geo %>% filter(is.na(continent) | is.na(country))
if (nrow(missing_geo)) {
  fwrite(missing_geo, file.path(OUT_DIR, "tables", "geography_missing_metadata_CHECK.csv"))
  stop(nrow(missing_geo),
       " model-score samples lack geographic metadata (see CHECK.csv).")
}
# Fail on conflicting sample-level labels: never silently use a mismatched CSV.
expected_geo_health <- case_when(
  dat_geo$geo_disease %in% c("healthy", "health", "control", "ctrl", "0") ~ "Healthy",
  dat_geo$geo_disease %in% c("disease", "diseased", "nonhealthy", "non-healthy", "case", "1") ~ "Disease",
  TRUE ~ NA_character_
)
# One sample can have no usable model outcome and is excluded later; do not
# turn that NA into an R `if (NA)` error. Check every KNOWN model outcome.
model_health <- as.character(dat_geo$health_status)
conflict_rows <- which(!is.na(model_health) &
                         (is.na(expected_geo_health) | expected_geo_health != model_health))
if (length(conflict_rows)) {
  data.table::fwrite(
    dat_geo[conflict_rows, c("sample_id", "disease", "geo_disease",
                             "dataset_split"), drop = FALSE],
    file.path(OUT_DIR, "tables", "metadata_health_status_conflicts_CHECK.csv")
  )
  stop(length(conflict_rows), " health-status mismatches between model scores",
       " and the matched metadata; inspect metadata_health_status_conflicts_CHECK.csv.")
}
# Exclude Oceania ONLY for geography panels and geography-specific tables.
# Keep the original dat/dat_use unchanged, so age/sex panels are unaffected.
geo_excluded_oceania <- sum(dat_geo$continent == "Oceania", na.rm = TRUE)
dat_geo <- dat_geo %>%
  filter(!is.na(OMWI), !is.na(health_status), continent %in% GEO_LEVELS) %>%
  mutate(continent = factor(continent, levels = GEO_LEVELS))
message("[geography] Oceania excluded from g/h/i (n=", geo_excluded_oceania,
        "); age/sex panels remain unchanged.")
message("[geography] Three-continent full-sample n=", nrow(dat_geo), "; countries=",
        n_distinct(dat_geo$country), "; studies=", n_distinct(dat_geo$geo_study_name))


# =========================================================
# 4. Statistical helper functions
# =========================================================
get_pairwise_wilcox <- function(df, group_var) {
  groups <- levels(droplevels(factor(df[[group_var]])))
  if (length(groups) < 2) return(data.frame())
  
  pairs <- combn(groups, 2, simplify = FALSE)
  out <- lapply(pairs, function(pair_now) {
    tmp <- df %>% filter(.data[[group_var]] %in% pair_now)
    tmp$group_tmp <- factor(tmp[[group_var]])
    p_now <- tryCatch(
      wilcox.test(OMWI ~ group_tmp, data = tmp, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    data.frame(
      group1 = pair_now[1],
      group2 = pair_now[2],
      p = p_now,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(out) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      p_label = p_to_star(p_adj)
    )
}

get_stratified_wilcox <- function(df, strata_var) {
  stat_df <- df %>%
    filter(!is.na(.data[[strata_var]]), !is.na(health_status), !is.na(OMWI)) %>%
    group_by(stratum = .data[[strata_var]]) %>%
    group_modify(function(.x, .y) {
      class_n <- table(.x$health_status)
      has_both <- all(c("Healthy", "Disease") %in% names(class_n))
      enough_n <- has_both && all(class_n[c("Healthy", "Disease")] >= MIN_PER_CLASS_FOR_VIOLIN)
      
      if (!enough_n) {
        return(tibble(
          n_healthy = ifelse("Healthy" %in% names(class_n), as.numeric(class_n["Healthy"]), 0),
          n_disease = ifelse("Disease" %in% names(class_n), as.numeric(class_n["Disease"]), 0),
          p = NA_real_
        ))
      }
      
      p_now <- tryCatch(
        wilcox.test(OMWI ~ health_status, data = .x, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
      
      tibble(
        n_healthy = as.numeric(class_n["Healthy"]),
        n_disease = as.numeric(class_n["Disease"]),
        p = p_now
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      p_label = p_to_star(p_adj)
    )
  
  stat_df
}

fit_interaction_lrt <- function(df, type = c("age", "sex"), adjusted = TRUE) {
  type <- match.arg(type)
  
  if (type == "age") {
    model_df <- df %>%
      filter(!is.na(disease01), !is.na(OMWI), !is.na(age_order)) %>%
      mutate(age_order = as.numeric(age_order))
    
    adjustment_terms <- character(0)
    if (adjusted && n_distinct(model_df$sex_group, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "sex_group")
    }
    if (adjusted && n_distinct(model_df$study_use, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "study_use")
    }
    
    reduced_terms <- c("OMWI", "age_order", adjustment_terms)
    full_terms <- c("OMWI * age_order", adjustment_terms)
  } else {
    model_df <- df %>%
      filter(!is.na(disease01), !is.na(OMWI), !is.na(sex_group))
    
    adjustment_terms <- character(0)
    if (adjusted && sum(!is.na(model_df$age_num)) >= 20) {
      adjustment_terms <- c(adjustment_terms, "age_num")
    }
    if (adjusted && n_distinct(model_df$study_use, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "study_use")
    }
    
    reduced_terms <- c("OMWI", "sex_group", adjustment_terms)
    full_terms <- c("OMWI * sex_group", adjustment_terms)
  }
  
  reduced_formula <- as.formula(paste("disease01 ~", paste(reduced_terms, collapse = " + ")))
  full_formula <- as.formula(paste("disease01 ~", paste(full_terms, collapse = " + ")))
  
  fit_reduced <- tryCatch(
    suppressWarnings(glm(reduced_formula, data = model_df, family = binomial())),
    error = function(e) NULL
  )
  fit_full <- tryCatch(
    suppressWarnings(glm(full_formula, data = model_df, family = binomial())),
    error = function(e) NULL
  )
  
  if (is.null(fit_reduced) || is.null(fit_full)) {
    return(tibble(
      analysis = type,
      adjusted = adjusted,
      n = nrow(model_df),
      p_interaction = NA_real_,
      model_reduced = paste(deparse(reduced_formula), collapse = ""),
      model_full = paste(deparse(full_formula), collapse = "")
    ))
  }
  
  lrt <- tryCatch(anova(fit_reduced, fit_full, test = "LRT"), error = function(e) NULL)
  p_int <- if (is.null(lrt)) NA_real_ else as.numeric(lrt$`Pr(>Chi)`[2])
  
  tibble(
    analysis = type,
    adjusted = adjusted,
    n = nrow(model_df),
    p_interaction = p_int,
    model_reduced = paste(deparse(reduced_formula), collapse = ""),
    model_full = paste(deparse(full_formula), collapse = "")
  )
}

calc_stratified_roc <- function(df, strata_var, strata_levels, colors_map,
                                min_per_class = MIN_PER_CLASS_FOR_AUC) {
  roc_objects <- list()
  roc_curve_list <- list()
  roc_stat_list <- list()
  
  for (stratum_now in strata_levels) {
    tmp <- df %>%
      filter(
        as.character(.data[[strata_var]]) == stratum_now,
        !is.na(disease01),
        !is.na(OMWI)
      )
    
    class_n <- table(tmp$disease01)
    n_healthy <- ifelse("0" %in% names(class_n), as.numeric(class_n["0"]), 0)
    n_disease <- ifelse("1" %in% names(class_n), as.numeric(class_n["1"]), 0)
    
    if (n_healthy < min_per_class || n_disease < min_per_class) {
      roc_stat_list[[stratum_now]] <- tibble(
        stratum = stratum_now,
        n_healthy = n_healthy,
        n_disease = n_disease,
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        sensitivity_at_OMWI0 = NA_real_,
        specificity_at_OMWI0 = NA_real_,
        balanced_accuracy_at_OMWI0 = NA_real_
      )
      next
    }
    
    roc_obj <- tryCatch(
      pROC::roc(
        response = tmp$disease01,
        predictor = -tmp$OMWI,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(roc_obj)) next
    
    set.seed(SET_SEED)
    auc_ci <- tryCatch(
      pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = ROC_BOOT_N),
      error = function(e) c(NA_real_, NA_real_, NA_real_)
    )
    
    roc_df <- tibble(
      specificity = as.numeric(roc_obj$specificities),
      sensitivity = as.numeric(roc_obj$sensitivities),
      stratum = stratum_now,
      fpr = 1 - specificity
    ) %>%
      arrange(fpr, sensitivity)
    
    roc_objects[[stratum_now]] <- roc_obj
    roc_curve_list[[stratum_now]] <- roc_df
    pred_disease_at_zero <- ifelse(tmp$OMWI < 0, 1, 0)
    sensitivity_at_zero <- mean(pred_disease_at_zero[tmp$disease01 == 1] == 1, na.rm = TRUE)
    specificity_at_zero <- mean(pred_disease_at_zero[tmp$disease01 == 0] == 0, na.rm = TRUE)
    balanced_accuracy_at_zero <- mean(c(sensitivity_at_zero, specificity_at_zero), na.rm = TRUE)
    
    roc_stat_list[[stratum_now]] <- tibble(
      stratum = stratum_now,
      n_healthy = n_healthy,
      n_disease = n_disease,
      auc = as.numeric(pROC::auc(roc_obj)),
      ci_low = as.numeric(auc_ci[1]),
      ci_high = as.numeric(auc_ci[3]),
      sensitivity_at_OMWI0 = sensitivity_at_zero,
      specificity_at_OMWI0 = specificity_at_zero,
      balanced_accuracy_at_OMWI0 = balanced_accuracy_at_zero
    )
  }
  
  stats_df <- bind_rows(roc_stat_list) %>%
    mutate(
      stratum = factor(stratum, levels = strata_levels),
      color = unname(colors_map[as.character(stratum)]),
      legend_label = ifelse(
        is.na(auc),
        paste0(as.character(stratum), ": insufficient n"),
        paste0(
          as.character(stratum), ": AUC ", sprintf("%.2f", auc),
          " (", sprintf("%.2f", ci_low), "–", sprintf("%.2f", ci_high), ")"
        )
      )
    )
  
  if (length(roc_curve_list) == 0) {
    curve_df <- tibble(
      specificity = numeric(0),
      sensitivity = numeric(0),
      stratum = factor(character(0), levels = strata_levels),
      fpr = numeric(0),
      legend_label = character(0)
    )
  } else {
    curve_df <- bind_rows(roc_curve_list) %>%
      mutate(
        stratum = factor(stratum, levels = strata_levels),
        legend_label = stats_df$legend_label[match(as.character(stratum), as.character(stats_df$stratum))]
      )
  }
  
  list(
    roc_objects = roc_objects,
    curve_df = curve_df,
    stats_df = stats_df
  )
}

compare_auc_delong <- function(roc_objects) {
  valid_names <- names(roc_objects)
  if (length(valid_names) < 2) return(tibble())
  
  pairs <- combn(valid_names, 2, simplify = FALSE)
  out <- lapply(pairs, function(pair_now) {
    test_now <- tryCatch(
      pROC::roc.test(
        roc_objects[[pair_now[1]]],
        roc_objects[[pair_now[2]]],
        method = "delong",
        paired = FALSE
      ),
      error = function(e) NULL
    )
    
    p_now <- if (is.null(test_now)) NA_real_ else as.numeric(test_now$p.value)
    
    tibble(
      group1 = pair_now[1],
      group2 = pair_now[2],
      p = p_now
    )
  })
  
  bind_rows(out) %>%
    mutate(p_adj = p.adjust(p, method = "BH"))
}

# =========================================================
# 5. Plot functions
# =========================================================
plot_overall_distribution <- function(df, group_var, group_levels, colors_map,
                                      x_title, p_label, tag_label,
                                      rotate_x = FALSE, draw_median_line = FALSE,
                                      pairwise_df = NULL, use_ggpubr = FALSE) {
  sub <- df %>%
    filter(!is.na(.data[[group_var]]), !is.na(OMWI)) %>%
    mutate(group_plot = factor(as.character(.data[[group_var]]), levels = group_levels)) %>%
    filter(!is.na(group_plot))
  
  if (n_distinct(sub$group_plot) < 2) return(NULL)
  
  median_df <- sub %>%
    group_by(group_plot) %>%
    summarise(median_OMWI = median(OMWI, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(sub, aes(x = group_plot, y = OMWI, fill = group_plot, color = group_plot)) +
    geom_violin(width = 0.88, trim = FALSE, linewidth = 0.65, alpha = 0.65) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.35, show.legend = FALSE) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.45
    ) +
    geom_point(
      data = median_df,
      aes(x = group_plot, y = median_OMWI),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.5
    ) +
    scale_fill_manual(values = colors_map, drop = FALSE) +
    scale_color_manual(values = colors_map, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = x_title, y = "OMWI", tag = tag_label) +
    theme_cell() +
    theme(
      legend.position = "none",
      axis.text.x = if (rotate_x) {
        element_text(angle = 25, hjust = 1, vjust = 1)
      } else {
        element_text(angle = 0, hjust = 0.5, vjust = 0.5)
      },
      plot.margin = margin(6, 7, 5.5, 5.5)
    )
  
  # Keep the original corner P labels for age and sex only.
  if (!is.null(p_label) && length(p_label) == 1L && !is.na(p_label)) {
    p <- p + annotate("text", x = Inf, y = Inf, label = p_label,
                      hjust = 1.04, vjust = 1.25, size = ANNOTATION_SIZE_MM)
  }
  
  # Geography and site: use precomputed BH-adjusted pairwise P values.
  # ggpubr draws brackets; it does not recompute unadjusted tests.
  # Non-significant comparisons are hidden, leaving stars only.
  # Keep the original p_to_star thresholds (maximum three stars).
  if (use_ggpubr) {
    if (!is.null(pairwise_df) && nrow(pairwise_df) > 0L) {
      pw <- pairwise_df %>%
        filter(is.finite(p_adj), p_adj < 0.05) %>%
        mutate(x1 = match(group1, group_levels),
               x2 = match(group2, group_levels)) %>%
        filter(!is.na(x1), !is.na(x2)) %>%
        arrange(x2 - x1, x1)
      if (nrow(pw) > 0L) {
        yr <- range(sub$OMWI, na.rm = TRUE)
        y_rng <- diff(yr)
        if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1
        pw$y.position <- yr[2] + (0.12 + 0.12 * (seq_len(nrow(pw)) - 1)) * y_rng
        pw$stars <- as.character(pw$p_label)
        p <- p + ggpubr::stat_pvalue_manual(
          pw, label = "stars", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.015,
          bracket.size = 0.4, size = ANNOTATION_SIZE_MM,
          inherit.aes = FALSE
        )
      }
    }
    return(p)
  }
  
  if (draw_median_line) {
    p <- p +
      geom_line(
        data = median_df,
        aes(x = group_plot, y = median_OMWI, group = 1),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.75
      )
  }
  
  if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
    pw_plot <- pairwise_df %>%
      filter(!is.na(p_adj), p_label != "ns", p_label != "") %>%
      mutate(
        x1 = match(group1, group_levels),
        x2 = match(group2, group_levels)
      ) %>%
      filter(!is.na(x1), !is.na(x2)) %>%
      arrange(x1, x2)
    
    if (nrow(pw_plot) > 0) {
      y_max <- max(sub$OMWI, na.rm = TRUE)
      y_min <- min(sub$OMWI, na.rm = TRUE)
      y_rng <- y_max - y_min
      if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1
      
      pw_plot <- pw_plot %>%
        mutate(
          y_bracket = y_max + (0.07 + 0.09 * (row_number() - 1)) * y_rng,
          y_label = y_bracket + 0.025 * y_rng,
          x_mid = (x1 + x2) / 2
        )
      
      p <- p +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x1, y = y_bracket, yend = y_bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x2, xend = x2, y = y_bracket, yend = y_bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_text(
          data = pw_plot,
          aes(x = x_mid, y = y_label, label = p_label),
          inherit.aes = FALSE,
          size = ANNOTATION_SIZE_MM,
          color = "black"
        )
    }
  }
  
  p
}

plot_stratified_violin <- function(df, strata_var, strata_levels, tag_label,
                                   include_single_class = FALSE) {
  sub <- df %>%
    filter(!is.na(.data[[strata_var]]), !is.na(health_status), !is.na(OMWI)) %>%
    mutate(
      stratum = factor(as.character(.data[[strata_var]]), levels = strata_levels),
      health_status = factor(health_status, levels = c("Healthy", "Disease"))
    ) %>%
    filter(!is.na(stratum))
  
  stat_df <- get_stratified_wilcox(sub, "stratum")
  keep_strata <- stat_df %>%
    filter(if (include_single_class) {
      n_healthy + n_disease >= MIN_PER_CLASS_FOR_VIOLIN
    } else {
      n_healthy >= MIN_PER_CLASS_FOR_VIOLIN & n_disease >= MIN_PER_CLASS_FOR_VIOLIN
    }) %>%
    pull(stratum) %>%
    as.character()
  
  sub <- sub %>%
    filter(as.character(stratum) %in% keep_strata) %>%
    mutate(stratum = droplevels(stratum))
  
  stat_df <- stat_df %>%
    filter(as.character(stratum) %in% keep_strata)
  
  if (nrow(sub) == 0) return(list(plot = NULL, stats = stat_df))
  
  y_stats <- sub %>%
    group_by(stratum) %>%
    summarise(
      y_max = max(OMWI, na.rm = TRUE),
      y_min = min(OMWI, na.rm = TRUE),
      y_rng = y_max - y_min,
      .groups = "drop"
    ) %>%
    mutate(y_rng = ifelse(y_rng <= 0, 1, y_rng))
  
  stat_plot <- stat_df %>%
    filter(n_healthy >= MIN_PER_CLASS_FOR_VIOLIN,
           n_disease >= MIN_PER_CLASS_FOR_VIOLIN) %>%
    left_join(y_stats, by = "stratum") %>%
    mutate(
      # Use numeric x positions so that the annotation layers do not
      # inherit x = health_status from the main ggplot data.
      x1 = 1,
      x2 = 2,
      x_mid = 1.5,
      y_bracket = y_max + 0.13 * y_rng,
      y_tip = y_bracket - 0.025 * y_rng,
      y_label = y_bracket + 0.025 * y_rng,
      p_label_plot = ifelse(p_label == "ns", "ns", p_label)
    )
  
  p <- ggplot(sub, aes(x = health_status, y = OMWI, fill = health_status, color = health_status)) +
    geom_violin(width = 0.88, trim = FALSE, linewidth = 0.65, alpha = 0.65) +
    geom_jitter(width = 0.06, size = 0.8, alpha = 0.30, show.legend = FALSE) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.45
    ) +
    # Draw the significance bracket manually. This avoids compatibility
    # problems in older ggpubr versions where stat_pvalue_manual may
    # inherit the parent health_status aesthetic.
    geom_segment(
      data = stat_plot,
      aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x1, xend = x1, y = y_bracket, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x2, xend = x2, y = y_bracket, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_text(
      data = stat_plot,
      aes(x = x_mid, y = y_label, label = p_label_plot),
      inherit.aes = FALSE,
      size = ANNOTATION_SIZE_MM,
      color = "black"
    ) +
    facet_grid(. ~ stratum, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = HD_COLORS, drop = FALSE) +
    scale_color_manual(values = HD_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = "OMWI", tag = tag_label) +
    theme_cell() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      panel.spacing.x = grid::unit(0.25, "lines"),
      plot.margin = margin(6, 6, 5.5, 5.5)
    )
  
  # Geography option: optionally retain a descriptive single-class facet.
  # Do not draw a significance bracket or pretend the AUC can be computed.
  if (include_single_class) {
    absent_class <- stat_df %>%
      filter(n_healthy == 0 | n_disease == 0) %>%
      left_join(y_stats, by = "stratum") %>%
      mutate(
        x_mid = 1.5,
        y_note = y_max + 0.14 * y_rng,
        note = ifelse(n_disease == 0, "No Disease samples", "No Healthy samples")
      )
    if (nrow(absent_class)) {
      p <- p + geom_text(
        data = absent_class,
        aes(x = x_mid, y = y_note, label = note),
        inherit.aes = FALSE,
        size = ANNOTATION_SIZE_MM * 0.85,
        color = "grey30"
      )
    }
  }
  
  list(plot = p, stats = stat_df)
}

plot_roc_panel <- function(roc_result, colors_map, interaction_label, tag_label) {
  curve_df <- roc_result$curve_df
  stats_df <- roc_result$stats_df %>% filter(!is.na(auc))
  
  if (nrow(curve_df) == 0) return(NULL)
  
  legend_values <- unname(colors_map[as.character(stats_df$stratum)])
  names(legend_values) <- stats_df$legend_label
  
  curve_df <- curve_df %>%
    mutate(legend_label = factor(legend_label, levels = stats_df$legend_label))
  
  p <- ggplot(curve_df, aes(x = fpr, y = sensitivity, color = legend_label)) +
    geom_step(linewidth = 1.0, direction = "vh") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey70", linewidth = 0.6) +
    scale_color_manual(values = legend_values, drop = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "1 − Specificity",
      y = "Sensitivity",
      color = NULL,
      tag = tag_label
    ) +
    theme_cell() +
    theme(
      legend.position = c(0.67, 0.27),
      legend.background = element_rect(fill = scales::alpha("white", 0.82), color = NA),
      legend.key.height = grid::unit(0.8, "lines"),
      legend.spacing.y = grid::unit(0.05, "lines"),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 5.5, 5.5)
    )
  
  if (SHOW_INTERACTION_ON_ROC && is.character(interaction_label) &&
      length(interaction_label) == 1L && !is.na(interaction_label)) {
    p <- p + annotate("text", x = 0.98, y = 0.04, label = interaction_label,
                      hjust = 1, vjust = 0, size = ANNOTATION_SIZE_MM * 0.9)
  }
  p
}

# =========================================================
# 6. Overall age and sex tests
# =========================================================
age_overall_df <- dat_use %>%
  filter(!is.na(age_group), !is.na(OMWI))

age_trend_test <- tryCatch(
  suppressWarnings(cor.test(as.numeric(age_overall_df$age_group), age_overall_df$OMWI, method = "spearman", exact = FALSE)),
  error = function(e) NULL
)
age_trend_p <- if (is.null(age_trend_test)) NA_real_ else age_trend_test$p.value

sex_overall_df <- dat_use %>%
  filter(!is.na(sex_group), !is.na(OMWI))

sex_wilcox <- tryCatch(
  wilcox.test(OMWI ~ sex_group, data = sex_overall_df, exact = FALSE),
  error = function(e) NULL
)
sex_p <- if (is.null(sex_wilcox)) NA_real_ else sex_wilcox$p.value
# =========================================================
# Overall tests shown in panels a and d
# =========================================================
overall_host_test_summary <- bind_rows(
  tibble(
    variable = "Age category",
    comparison = "18–44 vs 45–59 vs ≥60",
    test = "Spearman rank correlation using ordered age groups",
    n = nrow(age_overall_df),
    statistic_name = "Spearman rho",
    statistic = if (is.null(age_trend_test)) {
      NA_real_
    } else {
      as.numeric(age_trend_test$estimate)
    },
    p_value = age_trend_p,
    figure_label = format_p_label(age_trend_p, "Trend P")
  ),
  tibble(
    variable = "Gender",
    comparison = "Male vs Female",
    test = "Two-sided Wilcoxon rank-sum test",
    n = nrow(sex_overall_df),
    statistic_name = "Wilcoxon W",
    statistic = if (is.null(sex_wilcox)) {
      NA_real_
    } else {
      as.numeric(sex_wilcox$statistic)
    },
    p_value = sex_p,
    figure_label = format_p_label(sex_p, "P")
  )
)

print(overall_host_test_summary)
age_pairwise <- get_pairwise_wilcox(age_overall_df, "age_group")
sex_pairwise <- get_pairwise_wilcox(sex_overall_df, "sex_group")

# =========================================================
# 7. Stratified ROC analyses and heterogeneity tests
# =========================================================
age_levels <- c("18–44", "45–59", "≥60")
sex_levels <- c("Male", "Female")

age_roc <- calc_stratified_roc(
  df = dat_use,
  strata_var = "age_group",
  strata_levels = age_levels,
  colors_map = AGE_COLORS
)

sex_roc <- calc_stratified_roc(
  df = dat_use,
  strata_var = "sex_group",
  strata_levels = sex_levels,
  colors_map = SEX_COLORS
)

# Geography uses ALL data even when ANALYSIS_SCOPE="test" for age/sex.
# The geo ROC is a pooled, descriptive analysis, NOT independent validation.
geo_kw <- tryCatch(
  kruskal.test(OMWI ~ continent, data = dat_geo),
  error = function(e) NULL
)
geo_kw_p <- if (is.null(geo_kw)) NA_real_ else geo_kw$p.value
geo_pairwise <- get_pairwise_wilcox(dat_geo, "continent")
geo_roc <- calc_stratified_roc(
  df = dat_geo,
  strata_var = "continent",
  strata_levels = GEO_LEVELS,
  colors_map = GEO_COLORS
)
geo_delong <- compare_auc_delong(geo_roc$roc_objects) %>%
  mutate(analysis = "geography")

age_delong <- compare_auc_delong(age_roc$roc_objects) %>% mutate(analysis = "age")
sex_delong <- compare_auc_delong(sex_roc$roc_objects) %>% mutate(analysis = "sex")

interaction_tests <- bind_rows(
  fit_interaction_lrt(dat_use, type = "age", adjusted = FALSE),
  fit_interaction_lrt(dat_use, type = "age", adjusted = TRUE),
  fit_interaction_lrt(dat_use, type = "sex", adjusted = FALSE),
  fit_interaction_lrt(dat_use, type = "sex", adjusted = TRUE)
)

age_p_interaction <- interaction_tests %>%
  filter(analysis == "age", adjusted) %>%
  pull(p_interaction)
if (length(age_p_interaction) == 0 || is.na(age_p_interaction[1])) {
  age_p_interaction <- interaction_tests %>%
    filter(analysis == "age", !adjusted) %>%
    pull(p_interaction)
}
age_p_interaction <- age_p_interaction[1]

sex_p_interaction <- interaction_tests %>%
  filter(analysis == "sex", adjusted) %>%
  pull(p_interaction)
if (length(sex_p_interaction) == 0 || is.na(sex_p_interaction[1])) {
  sex_p_interaction <- interaction_tests %>%
    filter(analysis == "sex", !adjusted) %>%
    pull(p_interaction)
}
sex_p_interaction <- sex_p_interaction[1]

age_interaction_label <- format_p_label(age_p_interaction, "Interaction trend P")
sex_interaction_label <- format_p_label(sex_p_interaction, "Interaction P")

# =========================================================
# 7.0 Input inclusion audit (missing ages / sexes are NOT imputed)
# =========================================================
input_audit <- dat %>%
  group_by(dataset_split) %>%
  summarise(
    n_total = n(),
    n_valid_OMWI = sum(!is.na(OMWI)),
    n_known_health_status = sum(!is.na(health_status)),
    n_numeric_adult_age = sum(!is.na(age_group)),
    n_known_sex = sum(!is.na(sex_group)),
    .groups = "drop"
  )
print(input_audit)

# =========================================================
# 7.1 Subgroup sample-count tables
# Define these before figure rendering so they remain available even if
# a graphics device or plot layer fails later.
# =========================================================
age_counts <- dat_use %>%
  filter(!is.na(age_group), !is.na(health_status), !is.na(OMWI)) %>%
  count(age_group, health_status, name = "n") %>%
  complete(
    age_group = factor(age_levels, levels = age_levels, ordered = TRUE),
    health_status = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")),
    fill = list(n = 0)
  )

sex_counts <- dat_use %>%
  filter(!is.na(sex_group), !is.na(health_status), !is.na(OMWI)) %>%
  count(sex_group, health_status, name = "n") %>%
  complete(
    sex_group = factor(sex_levels, levels = sex_levels),
    health_status = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")),
    fill = list(n = 0)
  )

geo_counts <- dat_geo %>%
  count(continent, health_status, name = "n") %>%
  complete(
    continent = factor(GEO_LEVELS, levels = GEO_LEVELS),
    health_status = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")),
    fill = list(n = 0)
  )
geo_country_counts <- dat_geo %>%
  count(continent, country, dataset_split, health_status, name = "n")
geo_study_counts <- dat_geo %>%
  count(continent, country, geo_study_name, dataset_split, health_status, name = "n")
geo_coverage <- dat_geo %>%
  group_by(continent) %>%
  summarise(n_samples = n(), n_studies = n_distinct(geo_study_name),
            n_countries = n_distinct(country), .groups = "drop")

# =========================================================
# 8. Build panels
# =========================================================
p_age_overall <- plot_overall_distribution(
  df = dat_use,
  group_var = "age_group",
  group_levels = age_levels,
  colors_map = AGE_COLORS,
  x_title = "Age category",
  p_label = format_p_label(age_trend_p, "Trend P"),
  tag_label = "a",
  rotate_x = TRUE,
  draw_median_line = TRUE,
  pairwise_df = age_pairwise
)

age_stratified <- plot_stratified_violin(
  df = dat_use,
  strata_var = "age_group",
  strata_levels = age_levels,
  tag_label = "b"
)
p_age_stratified <- age_stratified$plot

p_age_roc <- plot_roc_panel(
  roc_result = age_roc,
  colors_map = AGE_COLORS,
  interaction_label = age_interaction_label,
  tag_label = "c"
)

p_sex_overall <- plot_overall_distribution(
  df = dat_use,
  group_var = "sex_group",
  group_levels = sex_levels,
  colors_map = SEX_COLORS,
  x_title = "Gender",
  p_label = format_p_label(sex_p, "P"),
  tag_label = "d",
  rotate_x = FALSE,
  draw_median_line = FALSE,
  pairwise_df = sex_pairwise
)

sex_stratified <- plot_stratified_violin(
  df = dat_use,
  strata_var = "sex_group",
  strata_levels = sex_levels,
  tag_label = "e"
)
p_sex_stratified <- sex_stratified$plot

p_sex_roc <- plot_roc_panel(
  roc_result = sex_roc,
  colors_map = SEX_COLORS,
  interaction_label = sex_interaction_label,
  tag_label = "f"
)

# Third row: original plot functions and drawing styles; distinct subgroup colors.
p_geo_overall <- plot_overall_distribution(
  df = dat_geo,
  group_var = "continent",
  group_levels = GEO_LEVELS,
  colors_map = GEO_COLORS,
  x_title = "Continent",
  p_label = NULL,
  tag_label = "g",
  rotate_x = TRUE,
  draw_median_line = FALSE,
  pairwise_df = geo_pairwise,
  use_ggpubr = TRUE
)

geo_stratified <- plot_stratified_violin(
  df = dat_geo,
  strata_var = "continent",
  strata_levels = GEO_LEVELS,
  tag_label = "h",
  include_single_class = FALSE
)
p_geo_stratified <- geo_stratified$plot

p_geo_roc <- plot_roc_panel(
  roc_result = geo_roc,
  colors_map = GEO_COLORS,
  interaction_label = NA_character_,
  tag_label = "i"
)

# Fourth row: oral sampling sites; same frozen model scores and metadata.
# Reuse the same model-ID-matched metadata subset as the geography panels.
# Metadata-only, unscored samples never enter the site plots.
# No second potentially conflicting site-lookup CSV is read.
site_normalize <- function(x) {
  x <- tolower(trimws(safe_text(x)))
  z <- rep(NA_character_, length(x))
  z[grepl("oral[ _-]*rinse|mouth[ _-]*rinse|mouthwash|oral wash|rinse", x)] <- "Oral rinse"
  z[grepl("pharyn|throat", x)] <- "Pharyngeal swab"
  z[grepl("tongue", x)] <- "Tongue dorsum"
  z[grepl("saliva|salivary", x)] <- "Saliva"
  z[grepl("dental|tooth|teeth|plaque|gingiv|periodontal|subging|supraging", x)] <- "Dental"
  z
}

if (!"body_site" %in% names(geo_meta)) {
  stop("Current metadata has no body_site; cannot generate oral-site panels.")
}
site_lookup <- geo_meta %>%
  transmute(sample_id = trimws(safe_text(sample_id)),
            site = site_normalize(body_site))
if (anyNA(site_lookup$site)) {
  stop(sum(is.na(site_lookup$site)),
       " model samples have an unrecognized oral sampling site. Check body_site.")
}
if (anyDuplicated(site_lookup$sample_id)) {
  stop("Duplicate sample IDs in current site lookup.")
}

# Use dat already created in the ORIGINAL script (from the frozen model).
# This preserves exactly the same sample scores, labels and split IDs used
# for the existing age/sex/geography figures.
site_dat_all <- dat %>%
  mutate(sample_id = trimws(safe_text(sample_id))) %>%
  left_join(site_lookup, by = "sample_id")
if (anyNA(site_dat_all$site)) {
  stop(sum(is.na(site_dat_all$site)), " model samples have no oral-site assignment in CURRENT metadata.",
       "\nFirst missing IDs: ",
       paste(head(site_dat_all$sample_id[is.na(site_dat_all$site)], 12), collapse = ", "))
}
if (any(!site_dat_all$site %in% c("Saliva", "Tongue dorsum", "Pharyngeal swab", "Oral rinse"))) {
  stop("Unexpected site, possibly dental, found among included model samples.")
}
if (anyDuplicated(site_dat_all$sample_id)) stop("Oral site join generated duplicate sample IDs")

SITE_LEVELS <- c("Saliva", "Tongue dorsum", "Pharyngeal swab")
SITE_COLORS <- c(
  "Saliva" = "#669B7A",        # EXACT sage green from AGE_COLORS
  "Tongue dorsum" = "#D7A650", # EXACT ochre from AGE_COLORS
  "Pharyngeal swab" = "#9879B5" # EXACT purple from AGE_COLORS
)
# Reuse EXACT SAME subsite palette in panels a and c (do not recolor).
SITE_ROC_COLORS <- SITE_COLORS

# Both a and b use pooled, descriptive data; rinse excluded from BOTH.
site_dat_abc <- site_dat_all %>%
  filter(site %in% SITE_LEVELS, !is.na(OMWI)) %>%
  mutate(site = factor(site, levels = SITE_LEVELS))
if (nrow(site_dat_abc) == 0L) stop("No saliva/tongue/pharyngeal samples")

# j: Overall OMWI distributions across the three subsites, without splitting
# by training/testing or by health status; replicate existing panels a/d/g.
site_kw <- stats::kruskal.test(OMWI ~ site, data = site_dat_abc)
site_pairwise <- get_pairwise_wilcox(site_dat_abc, "site")
p_site_a <- plot_overall_distribution(
  df = site_dat_abc,
  group_var = "site",
  group_levels = SITE_LEVELS,
  colors_map = SITE_COLORS,
  x_title = "Oral sampling site",
  p_label = NULL,
  tag_label = "j",
  rotate_x = TRUE,
  draw_median_line = FALSE,
  pairwise_df = site_pairwise,
  use_ggpubr = TRUE
)

# k: Healthy vs disease WITHIN each subsite, again pooled discovery+external;
# replicate the original b/e/h blue/red violin + boxes + P-value brackets.
# Pharyngeal swab has few healthy samples; the plot and sample count table
# retain that limitation. Original function requires >=5 per group for P.
site_b_result <- plot_stratified_violin(
  df = site_dat_abc,
  strata_var = "site",
  strata_levels = SITE_LEVELS,
  tag_label = "k",
  include_single_class = FALSE
)
p_site_b <- site_b_result$plot

# l: 3 ROC curves, one per oral subsite, each pooling TRAINING + EXTERNAL.
# Same data and same site grouping as panel b, without any extra row duplication.
# Original ROC convention: disease01=1 for non-healthy, predictor=-OMWI.
# Pharyngeal swabs have only 6 healthy samples in the pooled dataset; the
# general helper defaults to >=10/class, so use a site-specific >=2/class here.
# This does NOT affect any previously calculated age/sex/geography ROC curves.
site_roc_result <- calc_stratified_roc(
  df = site_dat_abc %>% filter(!is.na(disease01)),
  strata_var = "site",
  strata_levels = SITE_LEVELS,
  colors_map = SITE_ROC_COLORS,
  min_per_class = 2L
)
if (length(site_roc_result$roc_objects) != length(SITE_LEVELS)) {
  stop("Expected three pooled site-specific ROC curves. Check Healthy/Disease counts per site.")
}
p_site_c <- plot_roc_panel(
  roc_result = site_roc_result,
  colors_map = SITE_ROC_COLORS,
  interaction_label = NA_character_,
  tag_label = "l"
)
if (any(vapply(list(p_site_a, p_site_b, p_site_c), is.null, logical(1)))) {
  stop("One of the site j/k/l panels failed; inspect subgroup sample counts.")
}


plot_list <- list(
  p_age_overall,
  p_age_stratified,
  p_age_roc,
  p_sex_overall,
  p_sex_stratified,
  p_sex_roc,
  p_geo_overall,
  p_geo_stratified,
  p_geo_roc,
  p_site_a,
  p_site_b,
  p_site_c
)

if (any(vapply(plot_list, is.null, logical(1)))) {
  stop("At least one panel could not be generated. Check subgroup sample counts and metadata columns.")
}

# Figure assembly deliberately uses wrap_plots and does not use the patchwork ampersand operator.
figure_all <- patchwork::wrap_plots(
  plots = plot_list,
  ncol = 3,
  widths = c(1.05, 1.90, 1.12),
  heights = c(1, 1, 1, 1)
)
print(figure_all)
# =========================================================
# 9. Export figure and individual panels
# =========================================================
ggsave(
  filename = file.path(OUT_DIR, "Supplementary_Fig4_age_sex_geography_site_v5.pdf"),
  plot = figure_all,
  width = 22.0,
  height = 21.0,
  units = "in",
  device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
)

ggsave(
  filename = file.path(OUT_DIR, "Supplementary_Fig4_age_sex_geography_site_v5.png"),
  plot = figure_all,
  width = 22.0,
  height = 21.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

panel_names <- c(
  "a_age_overall",
  "b_age_stratified",
  "c_age_roc",
  "d_gender_overall",
  "e_gender_stratified",
  "f_gender_roc",
  "g_continent_overall",
  "h_continent_stratified",
  "i_continent_roc",
  "j_site_overall",
  "k_site_stratified",
  "l_site_roc"
)

panel_widths <- rep(c(6.0, 9.4, 6.3), 4L)
stopifnot(length(plot_list) == 12L, length(panel_names) == 12L,
          length(panel_widths) == 12L)

for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(OUT_DIR, "individual_panels", paste0(panel_names[i], ".pdf")),
    plot = plot_list[[i]],
    width = panel_widths[i],
    height = 5.0,
    units = "in",
    device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  )
  ggsave(
    filename = file.path(OUT_DIR, "individual_panels", paste0(panel_names[i], ".png")),
    plot = plot_list[[i]],
    width = panel_widths[i],
    height = 5.0,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# =========================================================
# 10. Export statistical tables
# =========================================================
fwrite(input_audit, file.path(OUT_DIR, "tables", "v3_input_inclusion_audit.csv"))
fwrite(age_counts, file.path(OUT_DIR, "tables", "age_group_sample_counts.csv"))
fwrite(sex_counts, file.path(OUT_DIR, "tables", "gender_sample_counts.csv"))
fwrite(geo_counts, file.path(OUT_DIR, "tables", "geography_continent_sample_counts.csv"))
fwrite(geo_country_counts, file.path(OUT_DIR, "tables", "geography_country_split_counts.csv"))
fwrite(geo_study_counts, file.path(OUT_DIR, "tables", "geography_study_split_counts.csv"))
fwrite(geo_coverage, file.path(OUT_DIR, "tables", "geography_coverage.csv"))
fwrite(geo_pairwise, file.path(OUT_DIR, "tables", "geography_overall_pairwise_wilcoxon.csv"))
fwrite(geo_stratified$stats, file.path(OUT_DIR, "tables", "geography_stratified_healthy_disease_wilcoxon.csv"))
fwrite(geo_roc$stats_df, file.path(OUT_DIR, "tables", "geography_fullsample_auc.csv"))
fwrite(geo_delong, file.path(OUT_DIR, "tables", "geography_auc_pairwise_delong.csv"))
fwrite(tibble(test = "Kruskal-Wallis", n = nrow(dat_geo),
              excluded_oceania_n = geo_excluded_oceania, p_value = geo_kw_p,
              shown_in_figure = FALSE),
       file.path(OUT_DIR, "tables", "geography_overall_continent_test.csv"))
fwrite(age_pairwise, file.path(OUT_DIR, "tables", "age_overall_pairwise_wilcoxon.csv"))
fwrite(sex_pairwise, file.path(OUT_DIR, "tables", "gender_overall_wilcoxon.csv"))
fwrite(age_stratified$stats, file.path(OUT_DIR, "tables", "age_stratified_healthy_disease_wilcoxon.csv"))
fwrite(sex_stratified$stats, file.path(OUT_DIR, "tables", "gender_stratified_healthy_disease_wilcoxon.csv"))
fwrite(age_roc$stats_df, file.path(OUT_DIR, "tables", "age_stratified_auc.csv"))
fwrite(sex_roc$stats_df, file.path(OUT_DIR, "tables", "gender_stratified_auc.csv"))
fwrite(age_delong, file.path(OUT_DIR, "tables", "age_auc_pairwise_delong.csv"))
fwrite(sex_delong, file.path(OUT_DIR, "tables", "gender_auc_pairwise_delong.csv"))
fwrite(interaction_tests, file.path(OUT_DIR, "tables", "OMWI_age_gender_interaction_tests.csv"))
fwrite(
  overall_host_test_summary,
  file.path(
    OUT_DIR,
    "tables",
    "overall_age_gender_tests_shown_in_figure.csv"
  )
)
cat("\n[DONE] Output directory: ", OUT_DIR, "\n", sep = "")
cat("Age discrimination test: ", age_interaction_label, "\n", sep = "")
cat("Sex discrimination test: ", sex_interaction_label, "\n", sep = "")
cat("Geography: Asia/Europe/Americas only, pooled in-sample training + external; descriptive AUC only.\n")
print(geo_roc$stats_df %>% select(stratum, n_healthy, n_disease, auc, ci_low, ci_high))

cat("Important: 'all' pools in-sample training scores and external data. ",
    "ANALYSIS_SCOPE='test' restricts age/sex only; geography/site remain pooled.\n", sep = "")



# =========================================================
# 11. Oral-site statistical audit (original split-specific tables)
# =========================================================
# ---------------------------------------------------------------------------
# Output complete audit: overall, discovery, external, and pooled descriptive.
# The three-site subgroup rows match ALL j/k/l panels. 'All sites' in the
# supplementary audit tables additionally includes oral rinse; it does NOT
# represent the population used in the site-specific three-curve ROC panel.
# ---------------------------------------------------------------------------
site_stat_one <- function(dd, cohort_label, site_label) {
  dd <- dd[is.finite(dd$OMWI) & !is.na(dd$disease01), , drop = FALSE]
  n_h <- sum(dd$disease01 == 0L)
  n_d <- sum(dd$disease01 == 1L)
  result <- data.frame(
    cohort = cohort_label, site = site_label, n = nrow(dd),
    n_healthy = n_h, n_nonhealthy = n_d,
    median_OMWI = if (nrow(dd)) median(dd$OMWI) else NA_real_,
    auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    wilcoxon_p = NA_real_, status = "Single outcome class",
    stringsAsFactors = FALSE
  )
  if (!n_h || !n_d) return(result)
  ro <- pROC::roc(response = dd$disease01, predictor = -dd$OMWI,
                  levels = c(0, 1), direction = "<", quiet = TRUE)
  result$auc <- as.numeric(pROC::auc(ro))
  if (min(n_h, n_d) >= 2L) {
    ci <- tryCatch(as.numeric(pROC::ci.auc(ro, method = "delong")),
                   error = function(e) c(NA_real_, NA_real_, NA_real_))
    result$ci_low <- ci[1]; result$ci_high <- ci[3]
  }
  result$wilcoxon_p <- tryCatch(stats::wilcox.test(
    dd$OMWI[dd$disease01 == 0L], dd$OMWI[dd$disease01 == 1L],
    exact = FALSE)$p.value, error = function(e) NA_real_)
  result$status <- "Both outcome classes"
  result
}
site_stat_table <- function(dd, cohort_label) {
  do.call(rbind, c(list(site_stat_one(dd, cohort_label, "All sites")),
                   lapply(c(SITE_LEVELS, "Oral rinse"), function(s) {
                     site_stat_one(dd[dd$site == s, , drop = FALSE], cohort_label, s)
                   })))
}
site_training_stats <- site_stat_table(
  site_dat_all[site_dat_all$dataset_split == "train", , drop = FALSE],
  "Discovery (in-sample)"
)
site_external_stats <- site_stat_table(
  site_dat_all[site_dat_all$dataset_split == "test", , drop = FALSE],
  "External validation"
)
site_pooled_stats <- site_stat_table(site_dat_all, "Overall (descriptive)")
site_all_stats <- dplyr::bind_rows(site_training_stats, site_external_stats, site_pooled_stats)
site_all_stats$estimate_type <- ifelse(
  site_all_stats$cohort == "Discovery (in-sample)", "Apparent (in-sample)",
  ifelse(site_all_stats$cohort == "External validation", "Independent external",
         "Pooled descriptive (includes in-sample scores)"))

site_write <- function(z, filename) {
  utils::write.csv(z, file.path(OUT_DIR, "tables", filename),
                   row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
site_write(site_pairwise, "site_overall_pairwise_wilcoxon.csv")
site_write(site_all_stats, "Table_all_discovery_external_pooled.csv")
site_write(site_training_stats, "Table_discovery.csv")
site_write(site_external_stats, "Table_external_validation.csv")
site_write(site_pooled_stats, "Table_pooled_descriptive.csv")
site_write(site_b_result$stats, "Table_site_within_group_wilcoxon.csv")
site_write(site_roc_result$stats_df, "Table_pooled_by_site_ROC_3_curves.csv")
site_write(
  site_dat_abc %>% count(site, health_status, dataset_split, name = "n"),
  "Table_three_site_sample_counts.csv"
)
site_write(data.frame(
  test = "Kruskal-Wallis, three sites pooled across splits",
  n = nrow(site_dat_abc), p_value = as.numeric(site_kw$p.value),
  shown_in_figure = FALSE),
  "Table_overall_three_site_Kruskal_Wallis.csv"
)

message("\n[Oral-site j/k/l DONE] Output: ", OUT_DIR)
message("Panel j/k: ONLY 3 sites, pooled discovery + external; NO oral rinse.")
message("Panel l: 3 pooled BY-SITE ROC curves; saliva/tongue/pharyngeal only; training + external within site.")
print(site_all_stats[, c("cohort", "site", "n", "n_healthy", "n_nonhealthy", "auc")], row.names = FALSE)

