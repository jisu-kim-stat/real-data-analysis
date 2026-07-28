#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- args_full[grepl("^--file=", args_full)]
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("work/real_data_analysis/R/00_download_data.R")
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
workspace_dir <- normalizePath(file.path(project_dir, "../.."), mustWork = TRUE)
cli <- commandArgs(trailingOnly = TRUE)
accept_meps <- "--accept-meps-terms" %in% cli

local_r_lib <- file.path(workspace_dir, "work", "Rlib")
dir.create(local_r_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_r_lib, .libPaths()))

install_if_missing <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install.packages(
      missing,
      lib = local_r_lib,
      repos = "https://cloud.r-project.org",
      dependencies = FALSE
    )
  }
}
install_if_missing(c("data.table", "foreign", "hdm", "gbm", "ggplot2"))
library(data.table)

super_dir <- file.path(
  workspace_dir, "work", "data", "raw", "superconductivity"
)
super_file <- file.path(super_dir, "train.csv")
if (!file.exists(super_file)) {
  dir.create(super_dir, recursive = TRUE, showWarnings = FALSE)
  zip_path <- file.path(super_dir, "superconductivity.zip")
  download.file(
    "https://archive.ics.uci.edu/static/public/464/superconductivty+data.zip",
    zip_path,
    mode = "wb"
  )
  unzip(zip_path, exdir = super_dir)
  unlink(zip_path)
}

meps_root <- file.path(
  workspace_dir, "work", "external", "cqr"
)
meps_raw_dir <- file.path(meps_root, "get_meps_data")
meps_dataset_dir <- file.path(meps_root, "datasets")
raw_csv <- file.path(meps_raw_dir, "h192.csv")
processed_csv <- file.path(meps_dataset_dir, "meps_21_reg.csv")
dir.create(meps_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meps_dataset_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_csv) || !file.exists(processed_csv)) {
  if (!accept_meps) {
    stop(
      paste(
        "MEPS HC-192 is a public-use file, but users must read and accept",
        "the AHRQ data-use terms before downloading it.",
        "Review https://meps.ahrq.gov/data_stats/download_data/pufs/h192/h192doc.shtml#DataA",
        "and rerun with --accept-meps-terms."
      )
    )
  }
}

if (!file.exists(raw_csv)) {
  zip_path <- file.path(meps_raw_dir, "h192ssp.zip")
  download.file(
    "https://meps.ahrq.gov/mepsweb/data_files/pufs/h192ssp.zip",
    zip_path,
    mode = "wb"
  )
  unzip(zip_path, exdir = meps_raw_dir)
  ssp_path <- file.path(meps_raw_dir, "h192.ssp")
  raw <- foreign::read.xport(ssp_path)
  fwrite(as.data.table(raw), raw_csv)
  unlink(c(zip_path, ssp_path))
}

if (!file.exists(processed_csv)) {
  d <- fread(raw_csv)
  d[, V1 := .I - 1L]
  d[, RACE := as.integer(HISPANX == 2 & RACEV2X == 1)]
  d <- d[PANEL == 21]

  rename_map <- c(
    FTSTU53X = "FTSTU", ACTDTY53 = "ACTDTY", HONRDC53 = "HONRDC",
    RTHLTH53 = "RTHLTH", MNHLTH53 = "MNHLTH", CHBRON53 = "CHBRON",
    JTPAIN53 = "JTPAIN", PREGNT53 = "PREGNT", WLKLIM53 = "WLKLIM",
    ACTLIM53 = "ACTLIM", SOCLIM53 = "SOCLIM", COGLIM53 = "COGLIM",
    EMPST53 = "EMPST", REGION53 = "REGION", MARRY53X = "MARRY",
    AGE53X = "AGE", POVCAT16 = "POVCAT", INSCOV16 = "INSCOV"
  )
  setnames(d, old = names(rename_map), new = unname(rename_map), skip_absent = TRUE)

  d <- d[
    REGION >= 0 & AGE >= 0 & MARRY >= 0 & ASTHDX >= 0
  ]
  categorical <- c(
    "REGION", "SEX", "MARRY", "FTSTU", "ACTDTY", "HONRDC", "RTHLTH",
    "MNHLTH", "HIBPDX", "CHDDX", "ANGIDX", "MIDX", "OHRTDX", "STRKDX",
    "EMPHDX", "CHBRON", "CHOLDX", "CANCERDX", "DIABDX", "JTPAIN",
    "ARTHDX", "ARTHTYPE", "ASTHDX", "ADHDADDX", "PREGNT", "WLKLIM",
    "ACTLIM", "SOCLIM", "COGLIM", "DFHEAR42", "DFSEE42", "ADSMOK42",
    "PHQ242", "EMPST", "POVCAT", "INSCOV"
  )
  valid_other <- setdiff(categorical, c("REGION", "MARRY", "ASTHDX"))
  d <- d[d[, rowSums(.SD < -1) == 0, .SDcols = valid_other]]
  use_cols <- c("OBTOTV16", "OPTOTV16", "ERTOT16", "IPNGTD16", "HHTOTD16")
  d <- d[d[, rowSums(.SD < 0) == 0, .SDcols = use_cols]]
  d[, UTILIZATION_reg := rowSums(.SD), .SDcols = use_cols]

  continuous <- c("AGE", "PCS42", "MCS42", "K6SUM42")
  keep <- c(
    "V1", continuous, "RACE", "UTILIZATION_reg", "PERWT16F", categorical
  )
  d <- d[, ..keep]
  d <- na.omit(d)

  one_hot <- lapply(categorical, function(v) {
    levels <- sort(unique(d[[v]]))
    out <- lapply(levels, function(level) as.integer(d[[v]] == level))
    names(out) <- paste0(v, "=", levels)
    as.data.table(out)
  })
  result <- cbind(
    d[, .(V1, AGE, RACE, PCS42, MCS42, K6SUM42, UTILIZATION_reg, PERWT16F)],
    rbindlist(list(do.call(cbind, one_hot)), use.names = TRUE)
  )
  fwrite(result, processed_csv)
}

cat("Data bootstrap complete.\n")
cat("Superconductivity:", super_file, "\n")
cat("MEPS raw:", raw_csv, "\n")
cat("MEPS processed:", processed_csv, "\n")
