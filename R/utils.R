#' Validate a subject count
#'
#' @param n Number of subjects.
#' @return `n` as an integer, invisibly.
#' @noRd
check_n <- function(n) {
  if (length(n) != 1L || is.na(n) || !is.numeric(n) || n < 1 ||
        n != floor(n)) {
    stop("`n` must be a single positive whole number.",
         call. = FALSE)
  }
  as.integer(n)
}

#' Validate an allocation probability
#'
#' @param p Probability of assignment to the deficient arm.
#' @param nm Argument name, for the error message.
#' @noRd
check_p <- function(p, nm = "p") {
  if (length(p) != 1L || is.na(p) || !is.numeric(p) ||
        p < 0 || p > 1) {
    stop(sprintf("`%s` must be a single probability in [0, 1].", nm),
         call. = FALSE)
  }
  as.numeric(p)
}

#' Validate covariate data supplied to a covariate-adaptive scheme
#'
#' Allocation procedures in the minimization family balance the
#' levels of discrete factors, so continuous columns are rejected
#' here rather than silently treated as having one level per
#' subject. Categorize before calling.
#'
#' @param covariates A data frame of factors or discrete columns,
#'   one row per subject in arrival order.
#' @return The covariates as a data frame of factors.
#' @noRd
check_covariates <- function(covariates) {
  if (!is.data.frame(covariates) || !ncol(covariates)) {
    stop("`covariates` must be a data frame with at least one ",
         "column.", call. = FALSE)
  }
  if (!nrow(covariates)) {
    stop("`covariates` must have at least one row.", call. = FALSE)
  }
  if (anyNA(covariates)) {
    stop("`covariates` must not contain missing values; a ",
         "sequential procedure cannot allocate a subject whose ",
         "balancing factors are unknown.", call. = FALSE)
  }
  out <- covariates
  for (j in seq_along(out)) {
    x <- out[[j]]
    if (is.numeric(x) && any(x != floor(x))) {
      stop("column '", names(out)[j], "' is continuous; ",
           "minimization balances factor levels, so continuous ",
           "covariates must be categorized first.", call. = FALSE)
    }
    out[[j]] <- factor(x)
  }
  out
}

#' Resolve per-covariate weights
#'
#' Weights express how much the design cares about balance in each
#' covariate. Unnamed weights are matched positionally; named
#' weights are matched by name and must cover every covariate.
#'
#' @param weights Numeric vector, or `NULL` for equal weights.
#' @param cov_names Character vector of covariate names.
#' @return A named numeric vector, one entry per covariate.
#' @noRd
resolve_weights <- function(weights, cov_names) {
  k <- length(cov_names)
  if (is.null(weights)) {
    return(stats::setNames(rep(1, k), cov_names))
  }
  if (!is.numeric(weights) || anyNA(weights)) {
    stop("`weights` must be a numeric vector without missing ",
         "values.", call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("`weights` must be non-negative.", call. = FALSE)
  }
  if (!is.null(names(weights))) {
    missing_nm <- setdiff(cov_names, names(weights))
    if (length(missing_nm)) {
      stop("`weights` is missing entries for: ",
           paste(missing_nm, collapse = ", "), call. = FALSE)
    }
    return(weights[cov_names])
  }
  if (length(weights) != k) {
    stop("`weights` has length ", length(weights), " but there ",
         "are ", k, " covariates.", call. = FALSE)
  }
  stats::setNames(weights, cov_names)
}

#' Imbalance between two arm counts
#'
#' @param n0,n1 Counts assigned to control and treatment.
#' @param measure One of `"diff"`, `"sd"`, or `"squared"`.
#' @noRd
imbalance_of <- function(n0, n1, measure = c("diff", "sd",
                                             "squared")) {
  measure <- match.arg(measure)
  switch(measure,
    diff = abs(n1 - n0),
    # The sample SD of the two counts is a monotone transform of
    # the absolute difference, and is the measure Pocock and
    # Simon used. It is retained for fidelity to the original.
    sd = stats::sd(c(n0, n1)),
    squared = (n1 - n0)^2
  )
}

#' Draw an assignment given a probability of the deficient arm
#'
#' Centralizes the tie rule: when the arms are equally balanced
#' every procedure here falls back to a fair coin, which is what
#' makes them all reduce to simple randomization when the
#' balancing criterion is uninformative.
#'
#' @param score1,score0 Imbalance that would result from assigning
#'   the arriving subject to treatment and to control.
#' @param p Probability of taking the arm that minimizes imbalance.
#' @noRd
biased_draw <- function(score1, score0, p) {
  if (score1 < score0) {
    prob1 <- p
  } else if (score1 > score0) {
    prob1 <- 1 - p
  } else {
    prob1 <- 0.5
  }
  as.integer(stats::runif(1) < prob1)
}
