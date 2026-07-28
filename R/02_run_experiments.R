#!/usr/bin/env Rscript

options(warn = 1)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- args_full[grepl("^--file=", args_full)]
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("work/real_data_analysis/R/02_run_experiments.R")
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(project_dir, "R", "tolerance_methods.R"))

required <- c("data.table", "gbm")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
library(data.table)

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

n_reps <- as.integer(get_arg("reps", "20"))
n_cores <- as.integer(get_arg("cores", "4"))
dataset_filter <- strsplit(get_arg("datasets", "all"), ",", fixed = TRUE)[[1]]
base_seed <- as.integer(get_arg("seed", "20260728"))

content <- 0.90
confidence_alpha <- 0.05
asr_tau <- 0.10
gbm_control <- list(
  n_trees = as.integer(get_arg("trees", "300")),
  interaction_depth = 3L,
  shrinkage = 0.04,
  bag_fraction = 0.75,
  min_obs_node = 20L
)

datasets <- readRDS(file.path(project_dir, "data", "processed", "datasets.rds"))
if (!identical(dataset_filter, "all")) {
  missing_ids <- setdiff(dataset_filter, names(datasets))
  if (length(missing_ids)) stop("Unknown datasets: ", paste(missing_ids, collapse = ", "))
  datasets <- datasets[dataset_filter]
}

result_dir <- file.path(project_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

run_one_split <- function(dataset, rep_id) {
  seed <- base_seed + match(dataset$id, names(datasets)) * 100000L + rep_id * 1000L
  set.seed(seed)
  n <- nrow(dataset$X)
  perm <- sample.int(n)
  n_train <- floor(n / 3)
  n_cal <- floor(n / 3)
  train_idx <- perm[seq_len(n_train)]
  cal_idx <- perm[n_train + seq_len(n_cal)]
  eval_idx <- perm[(n_train + n_cal + 1L):n]
  feature_keep <- apply(
    dataset$X[train_idx, , drop = FALSE],
    2,
    function(x) is.finite(var(x)) && var(x) > 0
  )
  X_train <- dataset$X[train_idx, feature_keep, drop = FALSE]
  X_cal <- dataset$X[cal_idx, feature_keep, drop = FALSE]
  X_eval <- dataset$X[eval_idx, feature_keep, drop = FALSE]

  models <- fit_nuisance_models(
    X_train = X_train,
    y_train = dataset$y_model[train_idx],
    content = content,
    asr_tau = asr_tau,
    seed = seed,
    gbm_control = gbm_control
  )
  pred_cal <- predict_nuisance(models, X_cal)
  pred_eval <- predict_nuisance(models, X_eval)
  built <- build_intervals(
    models = models,
    pred_cal = pred_cal,
    y_cal = dataset$y_model[cal_idx],
    pred_eval = pred_eval,
    content = content,
    confidence_alpha = confidence_alpha
  )

  y_eval_model <- dataset$y_model[eval_idx]
  y_eval_original <- dataset$y_original[eval_idx]
  scale_bin <- rank_bins(pred_eval$sigma, n_bins = 5L)

  method_rows <- list()
  group_rows <- list()
  score_rows <- list()

  for (method in names(built$intervals)) {
    int_model <- built$intervals[[method]]
    lower <- inverse_response(
      int_model[, "lower"], dataset$transform, dataset$lower_bound
    )
    upper <- inverse_response(
      int_model[, "upper"], dataset$transform, dataset$lower_bound
    )
    covered <- y_eval_original >= lower & y_eval_original <= upper
    lower_miss <- y_eval_original < lower
    upper_miss <- y_eval_original > upper
    width <- pmax(upper - lower, 0)

    method_rows[[method]] <- data.table(
      dataset = dataset$id,
      dataset_label = dataset$label,
      repetition = rep_id,
      method = method,
      n_train = length(train_idx),
      n_cal = length(cal_idx),
      n_eval = length(eval_idx),
      lambda = built$lambda,
      gamma = built$gamma,
      coverage = mean(covered),
      mean_width = mean(width),
      median_width = median(width),
      q90_width = unname(quantile(width, 0.90)),
      lower_miss = mean(lower_miss),
      upper_miss = mean(upper_miss),
      asr_a_minus = models$a_minus,
      asr_a_plus = models$a_plus
    )

    group_rows[[method]] <- rbindlist(lapply(1:5, function(g) {
      idx <- scale_bin == g
      data.table(
        dataset = dataset$id,
        repetition = rep_id,
        method = method,
        group_type = "predicted_scale_quintile",
        group = g,
        n = sum(idx),
        coverage = mean(covered[idx]),
        mean_width = mean(width[idx]),
        lower_miss = mean(lower_miss[idx]),
        upper_miss = mean(upper_miss[idx])
      )
    }))

    if (method %in% c("SR-TI", "ASR-TI", "CQR-TI")) {
      eval_score <- score_on_evaluation(
        method, models, pred_eval, y_eval_model
      )
      score_rows[[method]] <- rbindlist(lapply(1:5, function(g) {
        idx <- scale_bin == g
        data.table(
          dataset = dataset$id,
          repetition = rep_id,
          method = method,
          group = g,
          n = sum(idx),
          ks_to_marginal = ecdf_distance(eval_score[idx], eval_score)
        )
      }))
    }
  }

  list(
    metrics = rbindlist(method_rows),
    groups = rbindlist(group_rows),
    scores = rbindlist(score_rows)
  )
}

run_dataset <- function(dataset) {
  cat(
    sprintf(
      "[%s] n=%d p=%d reps=%d trees=%d\n",
      dataset$id, nrow(dataset$X), ncol(dataset$X), n_reps, gbm_control$n_trees
    )
  )
  runner <- function(r) {
    out <- run_one_split(dataset, r)
    cat(sprintf("[%s] completed repetition %d/%d\n", dataset$id, r, n_reps))
    out
  }
  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(seq_len(n_reps), runner, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_reps), runner)
  }
}

all_results <- lapply(datasets, run_dataset)
metrics <- rbindlist(lapply(all_results, function(x) rbindlist(lapply(x, `[[`, "metrics"))))
groups <- rbindlist(lapply(all_results, function(x) rbindlist(lapply(x, `[[`, "groups"))))
scores <- rbindlist(lapply(all_results, function(x) rbindlist(lapply(x, `[[`, "scores"))))

fwrite(metrics, file.path(result_dir, "split_metrics.csv"))
fwrite(groups, file.path(result_dir, "scale_group_metrics.csv"))
fwrite(scores, file.path(result_dir, "score_pivotality.csv"))
saveRDS(
  list(
    metrics = metrics,
    groups = groups,
    scores = scores,
    config = list(
      content = content,
      confidence_alpha = confidence_alpha,
      asr_tau = asr_tau,
      n_reps = n_reps,
      base_seed = base_seed,
      gbm_control = gbm_control
    )
  ),
  file.path(result_dir, "experiment_results.rds"),
  compress = "xz"
)

cat("Wrote experiment results to ", result_dir, "\n", sep = "")
