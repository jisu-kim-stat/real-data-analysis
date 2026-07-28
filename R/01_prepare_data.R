#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("work/real_data_analysis/R/01_prepare_data.R")
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
workspace_dir <- normalizePath(file.path(project_dir, "../.."), mustWork = TRUE)
processed_dir <- file.path(project_dir, "data", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

local_r_lib <- file.path(workspace_dir, "work", "Rlib")
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

required <- c("data.table", "hdm")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing R packages: ", paste(missing, collapse = ", "),
    ". Run the data-download/bootstrap step first."
  )
}

library(data.table)

make_cps <- function() {
  env <- new.env(parent = emptyenv())
  data("cps2012", package = "hdm", envir = env)
  cps <- get("cps2012", envir = env)

  rhs <- paste(
    c(
      "female", "widowed", "divorced", "separated", "nevermarried",
      "hsd08", "hsd911", "hsg", "cg", "ad", "mw", "so", "we",
      "exp1", "exp2"
    ),
    collapse = " + "
  )
  form <- as.formula(paste0("~ -1 + (", rhs, ")^2"))
  X <- model.matrix(form, data = cps)
  keep <- apply(X, 2, var) > 0
  X <- X[, keep, drop = FALSE]
  storage.mode(X) <- "double"

  education <- rep("other", nrow(cps))
  education[cps$hsd08 == 1 | cps$hsd911 == 1] <- "less_than_high_school"
  education[cps$hsg == 1] <- "high_school"
  education[cps$cg == 1] <- "college"
  education[cps$ad == 1] <- "advanced_degree"

  list(
    id = "cps_wage",
    label = "CPS 2012 hourly wage",
    X = X,
    y_model = as.numeric(cps$lnw),
    y_original = exp(as.numeric(cps$lnw)),
    transform = "log",
    lower_bound = 0,
    meta = data.frame(
      sex = ifelse(cps$female == 1, "female", "male"),
      education = education,
      experience = as.numeric(cps$exp1),
      survey_weight = as.numeric(cps$weight)
    ),
    source = "R package hdm::cps2012; DCP replication specification"
  )
}

make_superconductivity <- function() {
  path <- file.path(
    workspace_dir, "work", "data", "raw", "superconductivity", "train.csv"
  )
  if (!file.exists(path)) stop("Missing superconductivity file: ", path)
  d <- fread(path)
  y_name <- "critical_temp"
  feature_names <- setdiff(names(d), y_name)
  X <- as.matrix(d[, ..feature_names])
  storage.mode(X) <- "double"

  list(
    id = "superconductivity",
    label = "Superconductivity critical temperature",
    X = X,
    y_model = as.numeric(d[[y_name]]),
    y_original = as.numeric(d[[y_name]]),
    transform = "identity",
    lower_bound = 0,
    meta = data.frame(
      number_of_elements = as.numeric(d[["number_of_elements"]])
    ),
    source = "UCI Superconductivity Data"
  )
}

make_meps <- function() {
  processed_path <- file.path(
    workspace_dir, "work", "external", "cqr", "datasets", "meps_21_reg.csv"
  )
  raw_path <- file.path(
    workspace_dir, "work", "external", "cqr", "get_meps_data", "h192.csv"
  )
  if (!file.exists(processed_path)) stop("Missing CQR MEPS file: ", processed_path)
  if (!file.exists(raw_path)) stop("Missing raw MEPS file: ", raw_path)

  d <- fread(processed_path)
  raw <- fread(
    raw_path,
    select = c("PANEL", "TOTEXP16", "PERWT16F")
  )
  raw_row <- as.integer(d[["V1"]]) + 1L
  stopifnot(
    all(raw[["PANEL"]][raw_row] == 21),
    all(raw[["PERWT16F"]][raw_row] == d[["PERWT16F"]])
  )

  exclude <- c("V1", "UTILIZATION_reg", "PERWT16F")
  feature_names <- setdiff(names(d), exclude)
  X <- as.matrix(d[, ..feature_names])
  storage.mode(X) <- "double"

  health_cols <- paste0("RTHLTH=", 1:5)
  health_mat <- as.matrix(d[, ..health_cols])
  health_code <- max.col(health_mat, ties.method = "first")
  health_code[rowSums(health_mat) == 0] <- NA_integer_
  health_label <- c("excellent", "very_good", "good", "fair", "poor")

  common <- list(
    X = X,
    lower_bound = 0,
    meta = data.frame(
      age = as.numeric(d[["AGE"]]),
      race = ifelse(d[["RACE"]] == 1, "non_hispanic_white", "other"),
      sex = ifelse(d[["SEX=2"]] == 1, "female", "male"),
      self_rated_health = ifelse(
        is.na(health_code), "unknown", health_label[health_code]
      ),
      survey_weight = as.numeric(d[["PERWT16F"]])
    ),
    source = "AHRQ MEPS HC-192, Panel 21; CQR preprocessing cohort"
  )

  expenditure <- c(
    list(
      id = "meps_expenditure",
      label = "MEPS 2016 total medical expenditure",
      y_model = log1p(as.numeric(raw[["TOTEXP16"]][raw_row])),
      y_original = as.numeric(raw[["TOTEXP16"]][raw_row]),
      transform = "log1p"
    ),
    common
  )
  utilization <- c(
    list(
      id = "meps_utilization",
      label = "MEPS 2016 health-care utilization (sensitivity)",
      y_model = log1p(as.numeric(d[["UTILIZATION_reg"]])),
      y_original = as.numeric(d[["UTILIZATION_reg"]]),
      transform = "log1p"
    ),
    common
  )

  list(expenditure = expenditure, utilization = utilization)
}

datasets <- list(
  cps_wage = make_cps(),
  superconductivity = make_superconductivity()
)
meps <- make_meps()
datasets$meps_expenditure <- meps$expenditure
datasets$meps_utilization <- meps$utilization

descriptor <- rbindlist(lapply(datasets, function(d) {
  data.table(
    dataset = d$id,
    label = d$label,
    n = nrow(d$X),
    p = ncol(d$X),
    response_min = min(d$y_original),
    response_median = median(d$y_original),
    response_mean = mean(d$y_original),
    response_p90 = unname(quantile(d$y_original, 0.90)),
    response_p99 = unname(quantile(d$y_original, 0.99)),
    response_max = max(d$y_original),
    zero_fraction = mean(d$y_original == 0),
    transform = d$transform,
    source = d$source
  )
}))

saveRDS(datasets, file.path(processed_dir, "datasets.rds"), compress = "xz")
fwrite(descriptor, file.path(processed_dir, "dataset_descriptives.csv"))

cat("Prepared datasets:\n")
print(descriptor[, .(dataset, n, p, response_median, response_p99, zero_fraction)])
