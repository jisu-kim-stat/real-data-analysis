inverse_response <- function(x, transform, lower_bound = -Inf) {
  out <- switch(
    transform,
    identity = x,
    log = exp(x),
    log1p = expm1(x),
    stop("Unknown response transform: ", transform)
  )
  pmax(out, lower_bound)
}

empirical_quantile <- function(x, probability) {
  if (!length(x) || anyNA(x)) stop("Calibration scores must be finite and non-missing.")
  if (probability <= 0 || probability >= 1) {
    stop("Empirical quantile probability must lie strictly between 0 and 1.")
  }
  as.numeric(quantile(x, probs = probability, type = 1, names = FALSE))
}

rank_bins <- function(x, n_bins = 5L) {
  n <- length(x)
  pmin(n_bins, floor((rank(x, ties.method = "first") - 1) * n_bins / n) + 1L)
}

ecdf_distance <- function(x, reference) {
  grid <- sort(unique(c(x, reference)))
  fx <- findInterval(grid, sort(x)) / length(x)
  fr <- findInterval(grid, sort(reference)) / length(reference)
  max(abs(fx - fr))
}

fit_gbm <- function(
  X,
  y,
  distribution,
  seed,
  n_trees = 300L,
  interaction_depth = 3L,
  shrinkage = 0.04,
  bag_fraction = 0.75,
  min_obs_node = 20L
) {
  set.seed(seed)
  suppressWarnings(
    gbm::gbm.fit(
      x = X,
      y = y,
      distribution = distribution,
      n.trees = n_trees,
      interaction.depth = interaction_depth,
      shrinkage = shrinkage,
      bag.fraction = bag_fraction,
      n.minobsinnode = min_obs_node,
      nTrain = nrow(X),
      keep.data = FALSE,
      verbose = FALSE
    )
  )
}

predict_gbm <- function(model, X, n_trees) {
  as.numeric(predict(model, newdata = X, n.trees = n_trees))
}

normal_tolerance_factor <- function(content, confidence_alpha, df, mean_n) {
  z <- qnorm((1 + content) / 2)
  z * sqrt(df * (1 + 1 / mean_n) / qchisq(confidence_alpha, df = df))
}

fit_nuisance_models <- function(
  X_train,
  y_train,
  content,
  asr_tau,
  seed,
  gbm_control
) {
  n <- nrow(X_train)
  set.seed(seed + 11L)
  scale_idx <- sample.int(n, size = max(400L, floor(0.20 * n)))
  mean_idx <- setdiff(seq_len(n), scale_idx)

  mean_pre <- do.call(
    fit_gbm,
    c(
      list(
        X = X_train[mean_idx, , drop = FALSE],
        y = y_train[mean_idx],
        distribution = "gaussian",
        seed = seed + 101L
      ),
      gbm_control
    )
  )
  mu_scale <- predict_gbm(
    mean_pre, X_train[scale_idx, , drop = FALSE], gbm_control$n_trees
  )
  honest_residual <- y_train[scale_idx] - mu_scale

  residual_floor <- max(sd(y_train) * 1e-3, 1e-6)
  squared_residual <- honest_residual^2
  scale_cap <- unname(quantile(squared_residual, 0.995))
  scale_target <- pmin(squared_residual, scale_cap)
  scale_model <- do.call(
    fit_gbm,
    c(
      list(
        X = X_train[scale_idx, , drop = FALSE],
        y = scale_target,
        distribution = "gaussian",
        seed = seed + 202L,
        interaction_depth = max(2L, gbm_control$interaction_depth - 1L)
      ),
      gbm_control[setdiff(names(gbm_control), "interaction_depth")]
    )
  )

  mean_model <- do.call(
    fit_gbm,
    c(
      list(
        X = X_train,
        y = y_train,
        distribution = "gaussian",
        seed = seed + 303L
      ),
      gbm_control
    )
  )
  qlo_model <- do.call(
    fit_gbm,
    c(
      list(
        X = X_train,
        y = y_train,
        distribution = list(name = "quantile", alpha = (1 - content) / 2),
        seed = seed + 404L
      ),
      gbm_control
    )
  )
  qhi_model <- do.call(
    fit_gbm,
    c(
      list(
        X = X_train,
        y = y_train,
        distribution = list(name = "quantile", alpha = (1 + content) / 2),
        seed = seed + 505L
      ),
      gbm_control
    )
  )

  raw_sigma_scale <- sqrt(
    pmax(
      predict_gbm(
        scale_model, X_train[scale_idx, , drop = FALSE], gbm_control$n_trees
      ),
      residual_floor^2
    )
  )
  sigma_floor <- max(unname(quantile(raw_sigma_scale, 0.01)), residual_floor)
  raw_sigma_scale <- pmax(raw_sigma_scale, sigma_floor)
  scale_multiplier <- sqrt(mean((honest_residual / raw_sigma_scale)^2))
  z_honest <- honest_residual / (raw_sigma_scale * scale_multiplier)
  a_minus <- max(abs(empirical_quantile(z_honest, asr_tau)), 0.05)
  a_plus <- max(empirical_quantile(z_honest, 1 - asr_tau), 0.05)

  list(
    mean = mean_model,
    scale = scale_model,
    qlo = qlo_model,
    qhi = qhi_model,
    n_trees = gbm_control$n_trees,
    sigma_floor = sigma_floor,
    scale_multiplier = scale_multiplier,
    a_minus = a_minus,
    a_plus = a_plus,
    honest_residual_sd = sd(honest_residual),
    parametric_df = length(honest_residual) - 1L,
    parametric_mean_n = length(mean_idx)
  )
}

predict_nuisance <- function(models, X) {
  mu <- predict_gbm(models$mean, X, models$n_trees)
  sigma <- sqrt(
    pmax(
      predict_gbm(models$scale, X, models$n_trees),
      models$sigma_floor^2
    )
  )
  sigma <- pmax(sigma, models$sigma_floor) * models$scale_multiplier
  qlo <- predict_gbm(models$qlo, X, models$n_trees)
  qhi <- predict_gbm(models$qhi, X, models$n_trees)
  crossing <- qlo > qhi
  if (any(crossing)) {
    midpoint <- (qlo[crossing] + qhi[crossing]) / 2
    qlo[crossing] <- midpoint
    qhi[crossing] <- midpoint
  }
  list(mu = mu, sigma = sigma, qlo = qlo, qhi = qhi)
}

build_intervals <- function(
  models,
  pred_cal,
  y_cal,
  pred_eval,
  content,
  confidence_alpha
) {
  n_cal <- length(y_cal)
  lambda <- sqrt(log(2 / confidence_alpha) / (2 * n_cal))
  gamma <- content + lambda
  if (gamma >= 1) stop("C + lambda must be less than one.")

  z_cal <- (y_cal - pred_cal$mu) / pred_cal$sigma
  scores <- list(
    `SR-TI` = abs(z_cal),
    `ASR-TI` = pmax(
      -z_cal / models$a_minus,
      z_cal / models$a_plus
    ),
    `CQR-TI` = pmax(
      pred_cal$qlo - y_cal,
      y_cal - pred_cal$qhi
    )
  )
  thresholds <- vapply(scores, empirical_quantile, numeric(1), probability = gamma)

  intervals <- list(
    `SR-TI` = cbind(
      lower = pred_eval$mu - thresholds[["SR-TI"]] * pred_eval$sigma,
      upper = pred_eval$mu + thresholds[["SR-TI"]] * pred_eval$sigma
    ),
    `ASR-TI` = cbind(
      lower = pred_eval$mu -
        thresholds[["ASR-TI"]] * models$a_minus * pred_eval$sigma,
      upper = pred_eval$mu +
        thresholds[["ASR-TI"]] * models$a_plus * pred_eval$sigma
    ),
    `CQR-TI` = cbind(
      lower = pred_eval$qlo - thresholds[["CQR-TI"]],
      upper = pred_eval$qhi + thresholds[["CQR-TI"]]
    )
  )

  k <- normal_tolerance_factor(
    content = content,
    confidence_alpha = confidence_alpha,
    df = models$parametric_df,
    mean_n = models$parametric_mean_n
  )
  intervals[["Parametric-TI"]] <- cbind(
    lower = pred_eval$mu - k * models$honest_residual_sd,
    upper = pred_eval$mu + k * models$honest_residual_sd
  )

  list(
    intervals = intervals,
    thresholds = thresholds,
    lambda = lambda,
    gamma = gamma
  )
}

score_on_evaluation <- function(method, models, pred, y) {
  z <- (y - pred$mu) / pred$sigma
  switch(
    method,
    `SR-TI` = abs(z),
    `ASR-TI` = pmax(-z / models$a_minus, z / models$a_plus),
    `CQR-TI` = pmax(pred$qlo - y, y - pred$qhi),
    stop("Score diagnostic is not defined for ", method)
  )
}
