#' Covariate-adjusted response-adaptive allocation
#'
#' Uses accumulating outcome data as well as covariates: after a
#' burn-in allocated by fair coin, a working model is refitted as
#' subjects arrive, the arriving subject's success probability
#' under each arm is predicted from their covariates, and the
#' allocation probability is set by a target rule applied to those
#' predictions. The aim is to send subjects preferentially to the
#' arm that is working better *for people like them*, which is an
#' ethical objective rather than a balance one.
#'
#' This shifts the design's purpose relative to every other
#' procedure here. Those balance covariates and are agnostic to
#' outcomes; this one deliberately induces imbalance in response
#' to outcomes, and the resulting assignment sequence depends on
#' the responses, so standard tests require more care than under
#' covariate-adaptive randomization alone.
#'
#' The implementation assumes responses are available before the
#' next subject is randomized. That is a strong assumption and
#' frequently false in practice, since outcomes are often observed
#' long after enrolment; with delayed responses the procedure
#' behaves closer to its burn-in.
#'
#' @param covariates A data frame of covariates, one row per
#'   subject in arrival order.
#' @param response_fn A function of `(trt, covariate_row)`
#'   returning a single binary response, called once per subject
#'   with that subject's realized assignment.
#' @param target A function of `(p1, p0)`, the predicted success
#'   probabilities under treatment and control, returning the
#'   probability of assigning to treatment. The default is the
#'   RSIHR rule `sqrt(p1) / (sqrt(p1) + sqrt(p0))`, which
#'   minimizes expected failures for a fixed variance of the
#'   estimated treatment difference.
#' @param burn_in Subjects allocated by fair coin before adaptation
#'   begins. Must be large enough to fit the working model.
#'   Default 20.
#' @param clip Allocation probabilities are confined to
#'   `[clip, 1 - clip]`, so that no subject is assigned with
#'   near-certainty and every arm keeps accruing information.
#'   Default 0.1.
#' @return A list with `trt`, the integer assignment vector, and
#'   `response`, the realized responses.
#' @references Hu F, Rosenberger WF (2006). *The Theory of
#'   Response-Adaptive Randomization in Clinical Trials.* Wiley.
#'   Zhang L, Hu F (2009). A new family of covariate-adjusted
#'   response-adaptive designs and their properties. *Applied
#'   Mathematics, A Journal of Chinese Universities* 24:1-13.
#' @examples
#' set.seed(1)
#' cv <- data.frame(x = rnorm(80))
#' rf <- function(trt, row) {
#'   pr <- plogis(-0.2 + 0.8 * trt + 0.5 * row$x)
#'   rbinom(1, 1, pr)
#' }
#' out <- alloc_cara(cv, rf, burn_in = 30)
#' mean(out$trt)
#' @export
alloc_cara <- function(covariates, response_fn,
                       target = NULL, burn_in = 20L, clip = 0.1) {
  if (!is.data.frame(covariates) || !nrow(covariates)) {
    stop("`covariates` must be a data frame with at least one ",
         "row.", call. = FALSE)
  }
  if (!is.function(response_fn)) {
    stop("`response_fn` must be a function of (trt, row).",
         call. = FALSE)
  }
  clip <- check_p(clip, "clip")
  if (clip >= 0.5) {
    stop("`clip` must be below 0.5.", call. = FALSE)
  }
  if (is.null(target)) {
    target <- function(p1, p0) {
      s1 <- sqrt(max(p1, 0))
      s0 <- sqrt(max(p0, 0))
      if (s1 + s0 <= 0) 0.5 else s1 / (s1 + s0)
    }
  }
  n <- nrow(covariates)
  burn_in <- min(as.integer(burn_in), n)
  trt <- integer(n)
  resp <- numeric(n)
  cov_names <- names(covariates)
  fml <- stats::as.formula(paste(
    "resp ~ trt +", paste(cov_names, collapse = " + ")))

  for (i in seq_len(n)) {
    prob1 <- 0.5
    if (i > burn_in) {
      hist <- cbind(covariates[seq_len(i - 1L), , drop = FALSE],
                    trt = trt[seq_len(i - 1L)],
                    resp = resp[seq_len(i - 1L)])
      # Adapt only when both arms and both outcomes are present;
      # otherwise the fit is degenerate and a fair coin is the
      # honest fallback.
      ok <- length(unique(hist$trt)) == 2L &&
        length(unique(hist$resp)) == 2L
      if (ok) {
        fit <- tryCatch(
          suppressWarnings(stats::glm(fml, data = hist,
                                      family = stats::binomial())),
          error = function(e) NULL)
        if (!is.null(fit) && all(is.finite(stats::coef(fit)))) {
          nd <- covariates[rep(i, 2L), , drop = FALSE]
          nd$trt <- c(1L, 0L)
          pr <- tryCatch(
            stats::predict(fit, newdata = nd, type = "response"),
            error = function(e) NULL)
          if (!is.null(pr) && all(is.finite(pr))) {
            prob1 <- target(pr[[1L]], pr[[2L]])
          }
        }
      }
    }
    prob1 <- min(max(prob1, clip), 1 - clip)
    a <- as.integer(stats::runif(1) < prob1)
    trt[i] <- a
    r <- response_fn(a, covariates[i, , drop = FALSE])
    if (length(r) != 1L || is.na(r)) {
      stop("`response_fn` must return a single non-missing ",
           "value (subject ", i, ").", call. = FALSE)
    }
    resp[i] <- as.numeric(r)
  }
  list(trt = trt, response = resp)
}
