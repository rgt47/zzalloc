#' Pocock-Simon minimization
#'
#' For each arriving subject the procedure computes the imbalance
#' that would result from each assignment, summed over the levels
#' of the balancing factors that the subject actually has, and
#' assigns to the arm minimizing that total with probability `p`.
#' Balance is enforced on the *margins* rather than on the joint
#' cells, which is what lets minimization handle many factors
#' without the cell-count explosion that defeats stratification.
#'
#' Supplying `weights` gives the weighted variant, in which each
#' factor contributes to the score in proportion to how much the
#' design cares about balance in it. This is the natural
#' specification when the factors differ in prognostic strength,
#' since an imbalance matters in proportion to a factor's
#' association with the outcome, and it is the only member of this
#' family that can express a graded ranking: a stratified design
#' can include a factor or exclude it, but cannot weight it.
#'
#' @param covariates A data frame of discrete balancing factors,
#'   one row per subject in arrival order.
#' @param p Probability of assignment to the arm that minimizes
#'   imbalance. Default 0.8. `p = 0.5` recovers simple
#'   randomization; `p = 1` is deterministic minimization, which is
#'   Taves's original proposal and is the most predictable member
#'   of the family.
#' @param weights Optional per-covariate weights, named or
#'   positional. `NULL` gives equal weights.
#' @param measure Imbalance measure applied within each factor
#'   level: `"sd"` (Pocock and Simon's own), `"diff"`, or
#'   `"squared"`.
#' @return An integer vector with one assignment per row of
#'   `covariates`.
#' @references Pocock SJ, Simon R (1975). Sequential treatment
#'   assignment with balancing for prognostic factors in the
#'   controlled clinical trial. *Biometrics* 31(1):103-115.
#'   Taves DR (1974). Minimization: a new method of assigning
#'   patients to treatment and control groups. *Clinical
#'   Pharmacology and Therapeutics* 15(5):443-453.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = sample(c("F", "M"), 40, TRUE),
#'                  stage = sample(c("I", "II", "III"), 40, TRUE))
#' a <- alloc_pocock_simon(cv, p = 0.8)
#' # Weighted: balance stage four times as hard as sex.
#' w <- alloc_pocock_simon(cv, weights = c(sex = 1, stage = 4))
#' @export
alloc_pocock_simon <- function(covariates, p = 0.8, weights = NULL,
                               measure = c("sd", "diff",
                                           "squared")) {
  measure <- match.arg(measure)
  covariates <- check_covariates(covariates)
  p <- check_p(p)
  w <- resolve_weights(weights, names(covariates))
  n <- nrow(covariates)
  k <- ncol(covariates)

  # counts[[j]] is a 2-row matrix of assignments so far within
  # each level of covariate j: row 1 control, row 2 treatment.
  counts <- lapply(covariates, function(f)
    matrix(0L, nrow = 2L, ncol = nlevels(f),
           dimnames = list(c("0", "1"), levels(f))))
  lev <- vapply(covariates, as.integer, integer(n))
  if (n == 1L) lev <- matrix(lev, nrow = 1L)
  out <- integer(n)

  for (i in seq_len(n)) {
    s1 <- 0
    s0 <- 0
    for (j in seq_len(k)) {
      l <- lev[i, j]
      n0 <- counts[[j]][1L, l]
      n1 <- counts[[j]][2L, l]
      # Only the arriving subject's own levels enter the score;
      # levels they do not belong to are unaffected by where this
      # subject goes and would contribute the same constant to
      # both hypothetical totals.
      s1 <- s1 + w[[j]] * imbalance_of(n0, n1 + 1L, measure)
      s0 <- s0 + w[[j]] * imbalance_of(n0 + 1L, n1, measure)
    }
    a <- biased_draw(s1, s0, p)
    out[i] <- a
    for (j in seq_len(k)) {
      counts[[j]][a + 1L, lev[i, j]] <-
        counts[[j]][a + 1L, lev[i, j]] + 1L
    }
  }
  out
}

#' Hu-Hu joint balance minimization
#'
#' Extends Pocock-Simon by scoring three kinds of imbalance at
#' once: overall imbalance in the total sample size, marginal
#' imbalance within each factor level, and imbalance within the
#' subject's own joint stratum. The stratum term is what
#' distinguishes it, since a design balancing only margins can be
#' badly unbalanced within cells, which matters when the outcome
#' depends on a covariate interaction.
#'
#' Hu and Hu establish that the family attains bounded imbalance
#' on all three scales simultaneously under conditions on the
#' weights, which is not true of marginal balancing alone.
#'
#' @param covariates A data frame of discrete balancing factors.
#' @param p Probability of assignment to the arm minimizing the
#'   weighted score. Default 0.8.
#' @param weights Optional per-covariate weights for the marginal
#'   terms, named or positional.
#' @param overall_weight Weight on imbalance in the total sample
#'   size. Default 1.
#' @param stratum_weight Weight on imbalance within the subject's
#'   joint stratum. Default 1. Setting this to 0 and
#'   `overall_weight` to 0 recovers [alloc_pocock_simon()].
#' @param measure Imbalance measure. See [alloc_pocock_simon()].
#' @return An integer vector with one assignment per row.
#' @references Hu Y, Hu F (2012). Asymptotic properties of
#'   covariate-adaptive randomization. *Annals of Statistics*
#'   40(3):1794-1815.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = sample(c("F", "M"), 40, TRUE),
#'                  stage = sample(c("I", "II"), 40, TRUE))
#' alloc_hu_hu(cv, stratum_weight = 2)
#' @export
alloc_hu_hu <- function(covariates, p = 0.8, weights = NULL,
                        overall_weight = 1, stratum_weight = 1,
                        measure = c("sd", "diff", "squared")) {
  measure <- match.arg(measure)
  covariates <- check_covariates(covariates)
  p <- check_p(p)
  w <- resolve_weights(weights, names(covariates))
  for (nm in c("overall_weight", "stratum_weight")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop(sprintf("`%s` must be a single non-negative number.",
                   nm), call. = FALSE)
    }
  }
  n <- nrow(covariates)
  k <- ncol(covariates)

  counts <- lapply(covariates, function(f)
    matrix(0L, nrow = 2L, ncol = nlevels(f),
           dimnames = list(c("0", "1"), levels(f))))
  lev <- vapply(covariates, as.integer, integer(n))
  if (n == 1L) lev <- matrix(lev, nrow = 1L)
  stratum <- as.integer(interaction(covariates, drop = FALSE))
  ns <- max(stratum)
  scount <- matrix(0L, nrow = 2L, ncol = ns)
  tot <- c(0L, 0L)
  out <- integer(n)

  for (i in seq_len(n)) {
    s1 <- overall_weight * imbalance_of(tot[1L], tot[2L] + 1L,
                                        measure)
    s0 <- overall_weight * imbalance_of(tot[1L] + 1L, tot[2L],
                                        measure)
    for (j in seq_len(k)) {
      l <- lev[i, j]
      n0 <- counts[[j]][1L, l]
      n1 <- counts[[j]][2L, l]
      s1 <- s1 + w[[j]] * imbalance_of(n0, n1 + 1L, measure)
      s0 <- s0 + w[[j]] * imbalance_of(n0 + 1L, n1, measure)
    }
    st <- stratum[i]
    s1 <- s1 + stratum_weight *
      imbalance_of(scount[1L, st], scount[2L, st] + 1L, measure)
    s0 <- s0 + stratum_weight *
      imbalance_of(scount[1L, st] + 1L, scount[2L, st], measure)

    a <- biased_draw(s1, s0, p)
    out[i] <- a
    for (j in seq_len(k)) {
      counts[[j]][a + 1L, lev[i, j]] <-
        counts[[j]][a + 1L, lev[i, j]] + 1L
    }
    scount[a + 1L, st] <- scount[a + 1L, st] + 1L
    tot[a + 1L] <- tot[a + 1L] + 1L
  }
  out
}

#' Minimal sufficient balance
#'
#' Assigns by fair coin unless a covariate is *detectably*
#' unbalanced, and only then biases the coin toward the deficient
#' arm. Imbalance is judged by a test applied to each factor
#' level, and the coin is tilted only when a test falls below
#' `threshold`. Most allocations therefore remain purely random,
#' which preserves unpredictability while still preventing the
#' large imbalances that simple randomization occasionally
#' produces.
#'
#' This directly answers the predictability objection to
#' minimization: rather than tilting the coin at every allocation
#' to prevent imbalances that mostly would not occur, it
#' intervenes only when one has actually appeared.
#'
#' @param covariates A data frame of discrete balancing factors.
#' @param p Probability of the deficient arm when an imbalance is
#'   flagged. Default 0.8.
#' @param threshold Significance level below which a level is
#'   treated as imbalanced. Default 0.3, deliberately loose, since
#'   the aim is to detect imbalance early rather than to control
#'   an error rate.
#' @param burn_in Number of subjects allocated by fair coin before
#'   testing begins. Default 20.
#' @return An integer vector with one assignment per row.
#' @references Zhao W, Hill MD, Palesch Y (2015). Minimal
#'   sufficient balance: a new strategy to balance baseline
#'   covariates and preserve randomness of treatment allocation.
#'   *Statistical Methods in Medical Research* 24(6):989-1002.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = sample(c("F", "M"), 60, TRUE))
#' alloc_msb(cv)
#' @export
alloc_msb <- function(covariates, p = 0.8, threshold = 0.3,
                      burn_in = 20L) {
  covariates <- check_covariates(covariates)
  p <- check_p(p)
  threshold <- check_p(threshold, "threshold")
  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
        burn_in < 0) {
    stop("`burn_in` must be a single non-negative number.",
         call. = FALSE)
  }
  n <- nrow(covariates)
  k <- ncol(covariates)
  counts <- lapply(covariates, function(f)
    matrix(0L, nrow = 2L, ncol = nlevels(f),
           dimnames = list(c("0", "1"), levels(f))))
  lev <- vapply(covariates, as.integer, integer(n))
  if (n == 1L) lev <- matrix(lev, nrow = 1L)
  out <- integer(n)

  for (i in seq_len(n)) {
    prob1 <- 0.5
    if (i > burn_in) {
      # Vote across the subject's own levels: each flagged level
      # pushes toward its own deficient arm, and the direction
      # with more votes wins. A tie leaves the coin fair.
      votes <- 0L
      for (j in seq_len(k)) {
        l <- lev[i, j]
        n0 <- counts[[j]][1L, l]
        n1 <- counts[[j]][2L, l]
        tot <- n0 + n1
        if (tot < 2L) next
        # Two-sided binomial tail against an even split.
        pv <- stats::binom.test(min(n0, n1), tot, 0.5)$p.value
        if (pv < threshold) votes <- votes +
          if (n1 < n0) 1L else if (n1 > n0) -1L else 0L
      }
      if (votes > 0L) prob1 <- p
      if (votes < 0L) prob1 <- 1 - p
    }
    a <- as.integer(stats::runif(1) < prob1)
    out[i] <- a
    for (j in seq_len(k)) {
      counts[[j]][a + 1L, lev[i, j]] <-
        counts[[j]][a + 1L, lev[i, j]] + 1L
    }
  }
  out
}
