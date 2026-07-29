# Guo and Young (2024), "Approximate tolerance intervals for
# nonparametric regression models", Journal of Nonparametric Statistics.
#
# Exact implementation of the two-sided pointwise construction in:
#   - Equation (11): residual variance and Satterthwaite df
#   - Proposition 3.1 / Equation (14): pointwise k-factor
#
# The paper assumes the homoscedastic model
#   Y_i = f(x_i) + epsilon_i,  E(epsilon_i) = 0,
#   Var(epsilon_i) = sigma_epsilon^2,
# and a linear smoother f_hat(x) = ell_x^T Y.
#
# This wrapper uses R's GCV-selected cubic smoothing spline. It is therefore
# one-dimensional, matching the smoother used in the paper's simulations and
# data examples. P is content and gamma is confidence, so the paper's
# (0.90, 0.95) interval uses P = 0.90 and gamma = 0.95.

.guo2024_unpack_symmetric_band <- function(values, dimension, bandwidth = 4L) {
  stopifnot(
    length(values) == bandwidth * dimension,
    dimension >= 1L,
    bandwidth >= 1L
  )

  out <- matrix(0, nrow = dimension, ncol = dimension)
  for (offset in 0:(bandwidth - 1L)) {
    block <- values[offset * dimension + seq_len(dimension)]
    if (offset == 0L) {
      diag(out) <- block
    } else if (offset < dimension) {
      rows <- seq_len(dimension - offset)
      cols <- rows + offset
      out[cbind(rows, cols)] <- block[rows]
      out[cbind(cols, rows)] <- block[rows]
    }
  }
  out
}

.guo2024_spline_basis <- function(spline_fit, x_new) {
  fit <- spline_fit$fit
  x_new <- as.numeric(x_new)
  scaled_x <- (x_new - fit$min) / fit$range
  inside <- scaled_x >= 0 & scaled_x <= 1

  basis <- matrix(0, nrow = length(x_new), ncol = fit$nk)
  if (any(inside)) {
    basis[inside, ] <- splines::splineDesign(
      knots = fit$knot,
      x = scaled_x[inside],
      ord = 4L,
      derivs = 0L,
      outer.ok = TRUE
    )
  }

  left <- scaled_x < 0
  if (any(left)) {
    basis_left <- splines::splineDesign(
      knots = fit$knot,
      x = 0,
      ord = 4L,
      derivs = 0L,
      outer.ok = TRUE
    )
    slope_left <- splines::splineDesign(
      knots = fit$knot,
      x = 0,
      ord = 4L,
      derivs = 1L,
      outer.ok = TRUE
    )
    basis[left, ] <- matrix(
      basis_left,
      nrow = sum(left),
      ncol = fit$nk,
      byrow = TRUE
    ) + matrix(
      slope_left,
      nrow = sum(left),
      ncol = fit$nk,
      byrow = TRUE
    ) * scaled_x[left]
  }

  right <- scaled_x > 1
  if (any(right)) {
    basis_right <- splines::splineDesign(
      knots = fit$knot,
      x = 1,
      ord = 4L,
      derivs = 0L,
      outer.ok = TRUE
    )
    slope_right <- splines::splineDesign(
      knots = fit$knot,
      x = 1,
      ord = 4L,
      derivs = 1L,
      outer.ok = TRUE
    )
    delta <- scaled_x[right] - 1
    basis[right, ] <- matrix(
      basis_right,
      nrow = sum(right),
      ncol = fit$nk,
      byrow = TRUE
    ) + matrix(
      slope_right,
      nrow = sum(right),
      ncol = fit$nk,
      byrow = TRUE
    ) * delta
  }

  basis
}

.guo2024_trace <- function(x) {
  sum(diag(x))
}

guo2024_fit_smoothing_spline <- function(
  x,
  y,
  all_knots = FALSE,
  nknots = NULL,
  control_spar = list()
) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("x and y must have the same length.")
  if (length(x) < 4L) stop("At least four observations are required.")
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be finite.")
  }
  if (length(unique(x)) < 4L) {
    stop("At least four unique x values are required.")
  }

  fit_args <- list(
    x = x,
    y = y,
    cv = FALSE,
    all.knots = all_knots,
    keep.data = TRUE,
    keep.stuff = TRUE,
    control.spar = control_spar
  )
  if (!is.null(nknots)) fit_args$nknots <- nknots
  spline_fit <- do.call(stats::smooth.spline, fit_args)

  basis_train <- .guo2024_spline_basis(spline_fit, x)
  nk <- spline_fit$fit$nk
  cross_basis <- .guo2024_unpack_symmetric_band(
    spline_fit$auxM$XWX,
    dimension = nk
  )
  penalty <- .guo2024_unpack_symmetric_band(
    spline_fit$auxM$Sigma,
    dimension = nk
  )
  penalized_cross_basis <- cross_basis + spline_fit$lambda * penalty
  inverse_penalized_cross_basis <- solve(
    penalized_cross_basis,
    diag(nk)
  )

  beta <- drop(
    inverse_penalized_cross_basis %*% crossprod(basis_train, y)
  )
  fitted <- drop(basis_train %*% beta)
  reference_fitted <- stats::predict(spline_fit, x)$y
  fit_error <- max(abs(fitted - reference_fitted))
  if (!is.finite(fit_error) || fit_error > 1e-6) {
    stop(
      "Failed to reconstruct the smoothing-spline linear smoother; ",
      "maximum fitted-value discrepancy is ",
      signif(fit_error, 6),
      "."
    )
  }

  # L = B (B^T B + lambda Omega)^(-1) B^T.
  # The nonzero eigenvalues of L equal those of
  # H = (B^T B + lambda Omega)^(-1) B^T B.
  h_basis <- inverse_penalized_cross_basis %*% cross_basis
  h2 <- h_basis %*% h_basis
  h3 <- h2 %*% h_basis
  h4 <- h2 %*% h2
  tr_l <- .guo2024_trace(h_basis)
  tr_l2 <- .guo2024_trace(h2)
  tr_l3 <- .guo2024_trace(h3)
  tr_l4 <- .guo2024_trace(h4)

  # A = (I - L)^T (I - L).  The smoother matrix is symmetric here.
  trace_a <- length(y) - 2 * tr_l + tr_l2
  trace_a2 <- length(y) - 4 * tr_l + 6 * tr_l2 - 4 * tr_l3 + tr_l4
  if (trace_a <= 0 || trace_a2 <= 0) {
    stop("The reconstructed residual matrix has invalid traces.")
  }

  residual <- y - fitted
  sigma2_hat <- sum(residual^2) / trace_a
  nu <- trace_a^2 / trace_a2

  list(
    spline_fit = spline_fit,
    basis_train = basis_train,
    cross_basis = cross_basis,
    penalty = penalty,
    inverse_penalized_cross_basis = inverse_penalized_cross_basis,
    beta = beta,
    fitted = fitted,
    residual = residual,
    sigma2_hat = sigma2_hat,
    sigma_hat = sqrt(sigma2_hat),
    nu = nu,
    trace_a = trace_a,
    trace_a2 = trace_a2,
    lambda = spline_fit$lambda,
    spar = spline_fit$spar,
    effective_df = tr_l,
    fit_error = fit_error
  )
}

guo2024_equation14_probability <- function(
  k,
  ell_norm,
  nu,
  P,
  rel_tol = 1e-8,
  subdivisions = 300L
) {
  if (!is.finite(k) || k <= 0) stop("k must be positive and finite.")
  if (!is.finite(ell_norm) || ell_norm < 0) {
    stop("ell_norm must be nonnegative and finite.")
  }
  if (!is.finite(nu) || nu <= 0) stop("nu must be positive and finite.")
  if (!is.finite(P) || P <= 0 || P >= 1) stop("P must lie in (0, 1).")

  if (ell_norm <= sqrt(.Machine$double.eps)) {
    threshold <- nu * stats::qchisq(P, df = 1) / k^2
    return(stats::pchisq(threshold, df = nu, lower.tail = FALSE))
  }

  integrand <- function(z) {
    noncentral_quantile <- stats::qchisq(
      P,
      df = 1,
      ncp = (ell_norm * z)^2
    )
    threshold <- nu * noncentral_quantile / k^2
    sqrt(2 / pi) *
      stats::pchisq(threshold, df = nu, lower.tail = FALSE) *
      exp(-z^2 / 2)
  }

  stats::integrate(
    integrand,
    lower = 0,
    # The omitted half-normal probability above 12 is below 4e-33.
    # A finite upper bound also prevents numerical failures in qchisq()
    # caused by irrelevant quadrature evaluations at enormous z values.
    upper = 12,
    rel.tol = rel_tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )$value
}

guo2024_two_sided_k <- function(
  ell_norm,
  nu,
  P = 0.90,
  gamma = 0.95,
  rel_tol = 1e-8,
  root_tol = 1e-8,
  subdivisions = 300L
) {
  if (!is.finite(gamma) || gamma <= 0 || gamma >= 1) {
    stop("gamma must be the confidence level and must lie in (0, 1).")
  }

  if (ell_norm <= sqrt(.Machine$double.eps)) {
    return(
      sqrt(
        nu * stats::qchisq(P, df = 1) /
          stats::qchisq(1 - gamma, df = nu)
      )
    )
  }

  objective <- function(k) {
    guo2024_equation14_probability(
      k = k,
      ell_norm = ell_norm,
      nu = nu,
      P = P,
      rel_tol = rel_tol,
      subdivisions = subdivisions
    ) - gamma
  }

  lower <- 1e-6
  upper <- 2
  while (objective(upper) < 0) {
    upper <- upper * 2
    if (upper > 1e4) stop("Could not bracket the Equation (14) k-factor.")
  }

  stats::uniroot(
    objective,
    interval = c(lower, upper),
    tol = root_tol
  )$root
}

guo2024_predict_smoother <- function(model, x_new) {
  basis_new <- .guo2024_spline_basis(model$spline_fit, x_new)
  fitted <- drop(basis_new %*% model$beta)

  # ell_x = b(x)^T (B^T B + lambda Omega)^(-1) B^T.
  # Hence ||ell_x||^2 = b(x)^T K^(-1) B^T B K^(-1) b(x).
  ell_quadratic <- model$inverse_penalized_cross_basis %*%
    model$cross_basis %*%
    model$inverse_penalized_cross_basis
  ell_norm2 <- rowSums((basis_new %*% ell_quadratic) * basis_new)
  ell_norm <- sqrt(pmax(ell_norm2, 0))

  list(
    fitted = fitted,
    ell_norm = ell_norm,
    basis = basis_new
  )
}

guo2024_pointwise_ti <- function(
  x,
  y,
  x_new = x,
  P = 0.90,
  gamma = 0.95,
  all_knots = FALSE,
  nknots = NULL,
  control_spar = list(),
  k_cache_digits = 10L,
  rel_tol = 1e-8,
  root_tol = 1e-8,
  subdivisions = 300L
) {
  model <- guo2024_fit_smoothing_spline(
    x = x,
    y = y,
    all_knots = all_knots,
    nknots = nknots,
    control_spar = control_spar
  )
  prediction <- guo2024_predict_smoother(model, x_new)

  cache_key <- formatC(
    prediction$ell_norm,
    digits = k_cache_digits,
    format = "fg",
    flag = "#"
  )
  unique_key <- unique(cache_key)
  representative_norm <- vapply(
    unique_key,
    function(key) prediction$ell_norm[match(key, cache_key)],
    numeric(1)
  )
  unique_k <- vapply(
    representative_norm,
    guo2024_two_sided_k,
    numeric(1),
    nu = model$nu,
    P = P,
    gamma = gamma,
    rel_tol = rel_tol,
    root_tol = root_tol,
    subdivisions = subdivisions
  )
  k <- unname(unique_k[match(cache_key, unique_key)])

  half_width <- k * model$sigma_hat
  interval <- cbind(
    lower = prediction$fitted - half_width,
    upper = prediction$fitted + half_width
  )

  list(
    interval = interval,
    fitted = prediction$fitted,
    ell_norm = prediction$ell_norm,
    k = k,
    sigma_hat = model$sigma_hat,
    sigma2_hat = model$sigma2_hat,
    nu = model$nu,
    P = P,
    gamma = gamma,
    model = model
  )
}
