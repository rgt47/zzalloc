#' Imbalance diagnostics for a realized allocation
#'
#' Reports imbalance on the three scales the procedures in this
#' package target, so that a design can be judged on the scale it
#' was built for rather than on whichever is convenient.
#'
#' @param trt Integer vector of assignments, `0` or `1`.
#' @param covariates Optional data frame of balancing factors, one
#'   row per assignment.
#' @return A list with `overall`, the difference in arm sizes;
#'   `marginal`, a named list of per-level differences for each
#'   covariate; `max_marginal`, the largest absolute marginal
#'   difference; and `stratum`, per-cell differences for the joint
#'   cross-classification. The covariate entries are `NULL` when
#'   no covariates are supplied.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = sample(c("F", "M"), 40, TRUE))
#' a <- alloc_pocock_simon(cv)
#' allocation_imbalance(a, cv)
#' @export
allocation_imbalance <- function(trt, covariates = NULL) {
  trt <- check_trt(trt)
  overall <- sum(trt == 1L) - sum(trt == 0L)
  if (is.null(covariates)) {
    return(list(overall = overall, marginal = NULL,
                max_marginal = NA_real_, stratum = NULL))
  }
  covariates <- check_covariates(covariates)
  if (nrow(covariates) != length(trt)) {
    stop("`covariates` has ", nrow(covariates), " rows but `trt` ",
         "has ", length(trt), " entries.", call. = FALSE)
  }
  marg <- lapply(covariates, function(f) {
    tb <- table(factor(trt, levels = c(0L, 1L)), f)
    stats::setNames(as.integer(tb[2L, ] - tb[1L, ]), colnames(tb))
  })
  st <- interaction(covariates, drop = TRUE, sep = "|")
  tbs <- table(factor(trt, levels = c(0L, 1L)), st)
  strat <- stats::setNames(as.integer(tbs[2L, ] - tbs[1L, ]),
                           colnames(tbs))
  list(overall = overall, marginal = marg,
       max_marginal = max(abs(unlist(marg))), stratum = strat)
}

#' Response-weighted imbalance
#'
#' The difference between arms in the fitted linear predictor,
#' `sum_j beta_j * (xbar_1j - xbar_0j)`. Each covariate's
#' imbalance is weighted by its association with the outcome, so
#' an imbalance in a strong predictor counts for more than the
#' same imbalance in a weak one.
#'
#' This is the scale on which imbalance actually matters. A design
#' can look well balanced on covariate counts while being
#' materially unbalanced on this measure, if the imbalances happen
#' to fall in the covariates that predict the outcome, and the
#' converse also holds.
#'
#' @param trt Integer vector of assignments.
#' @param covariates Data frame of covariates, numeric or factors
#'   with two levels.
#' @param beta Named numeric vector of covariate effects on the
#'   outcome, one per column of `covariates`.
#' @return A single number: the response-weighted imbalance.
#'   Positive values mean the treatment arm has the higher fitted
#'   linear predictor.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = sample(0:1, 40, TRUE),
#'                  stage = sample(0:1, 40, TRUE))
#' a <- alloc_simple(40)
#' weighted_imbalance(a, cv, beta = c(sex = 0.2, stage = 1.1))
#' @export
weighted_imbalance <- function(trt, covariates, beta) {
  trt <- check_trt(trt)
  if (!is.data.frame(covariates)) {
    stop("`covariates` must be a data frame.", call. = FALSE)
  }
  if (nrow(covariates) != length(trt)) {
    stop("`covariates` has ", nrow(covariates), " rows but `trt` ",
         "has ", length(trt), " entries.", call. = FALSE)
  }
  b <- resolve_weights_signed(beta, names(covariates))
  if (!any(trt == 1L) || !any(trt == 0L)) {
    stop("both arms must be represented to compute an imbalance.",
         call. = FALSE)
  }
  total <- 0
  for (nm in names(covariates)) {
    x <- covariates[[nm]]
    if (is.factor(x)) {
      if (nlevels(x) != 2L) {
        stop("factor '", nm, "' has ", nlevels(x), " levels; the ",
             "response-weighted measure needs numeric or ",
             "two-level covariates.", call. = FALSE)
      }
      x <- as.integer(x) - 1L
    }
    if (!is.numeric(x)) {
      stop("column '", nm, "' is not numeric.", call. = FALSE)
    }
    total <- total + b[[nm]] *
      (mean(x[trt == 1L]) - mean(x[trt == 0L]))
  }
  total
}

#' @noRd
check_trt <- function(trt) {
  if (!length(trt) || anyNA(trt) ||
        !all(trt %in% c(0L, 1L, 0, 1))) {
    stop("`trt` must be a non-empty vector of 0s and 1s without ",
         "missing values.", call. = FALSE)
  }
  as.integer(trt)
}

#' Like resolve_weights but permits negative effects
#' @noRd
resolve_weights_signed <- function(beta, cov_names) {
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be a numeric vector without missing ",
         "values.", call. = FALSE)
  }
  if (!is.null(names(beta))) {
    missing_nm <- setdiff(cov_names, names(beta))
    if (length(missing_nm)) {
      stop("`beta` is missing entries for: ",
           paste(missing_nm, collapse = ", "), call. = FALSE)
    }
    return(beta[cov_names])
  }
  if (length(beta) != length(cov_names)) {
    stop("`beta` has length ", length(beta), " but there are ",
         length(cov_names), " covariates.", call. = FALSE)
  }
  stats::setNames(beta, cov_names)
}
