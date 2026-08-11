#' Efron's biased coin design
#'
#' When the arms are unequal the next subject goes to the smaller
#' arm with probability `p`; when they are equal a fair coin is
#' used. Efron (1971) introduced this as the first procedure to
#' trade a little unpredictability for a bound on imbalance, and it
#' remains the reference against which restricted procedures are
#' compared.
#'
#' @param n Number of subjects.
#' @param p Probability of assignment to the deficient arm.
#'   Default `2/3`, Efron's own choice. `p = 0.5` recovers simple
#'   randomization; `p = 1` gives deterministic alternation once
#'   the arms differ.
#' @return An integer vector of length `n`.
#' @references Efron B (1971). Forcing a sequential experiment to
#'   be balanced. *Biometrika* 58(3):403-417.
#' @examples
#' set.seed(1)
#' alloc_efron(10)
#' @export
alloc_efron <- function(n, p = 2 / 3) {
  n <- check_n(n)
  p <- check_p(p)
  out <- integer(n)
  n1 <- 0L
  n0 <- 0L
  for (i in seq_len(n)) {
    prob1 <- if (n1 < n0) p else if (n1 > n0) 1 - p else 0.5
    a <- as.integer(stats::runif(1) < prob1)
    out[i] <- a
    if (a == 1L) n1 <- n1 + 1L else n0 <- n0 + 1L
  }
  out
}

#' Wei's urn design
#'
#' An urn starts with `alpha` balls of each colour. For each
#' subject a ball is drawn and replaced, the corresponding arm is
#' assigned, and `beta` balls of the *opposite* colour are added.
#' The pull toward balance is therefore strong early, when the
#' added balls are a large fraction of the urn, and weakens as the
#' trial grows, so the procedure becomes progressively less
#' predictable. This is the opposite of the permuted block, whose
#' predictability is constant.
#'
#' @param n Number of subjects.
#' @param alpha Initial number of balls of each colour.
#' @param beta Balls of the opposite colour added after each draw.
#'   `beta = 0` gives simple randomization.
#' @return An integer vector of length `n`.
#' @references Wei LJ (1977). A class of designs for sequential
#'   clinical trials. *JASA* 72(358):382-386.
#' @examples
#' set.seed(1)
#' alloc_wei_urn(20, alpha = 1, beta = 1)
#' @export
alloc_wei_urn <- function(n, alpha = 1, beta = 1) {
  n <- check_n(n)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(beta) || length(beta) != 1L || beta < 0) {
    stop("`beta` must be a single non-negative number.",
         call. = FALSE)
  }
  out <- integer(n)
  n1 <- 0L
  n0 <- 0L
  for (i in seq_len(n)) {
    # Balls favouring arm 1 are the initial alpha plus beta for
    # every subject already sent to arm 0.
    w1 <- alpha + beta * n0
    w0 <- alpha + beta * n1
    a <- as.integer(stats::runif(1) < w1 / (w1 + w0))
    out[i] <- a
    if (a == 1L) n1 <- n1 + 1L else n0 <- n0 + 1L
  }
  out
}

#' Smith's generalized biased coin design
#'
#' The probability of the deficient arm varies smoothly with how
#' unbalanced the trial is, rather than jumping as in Efron's
#' design. With counts `n0` and `n1` the probability of assignment
#' to treatment is `n0^rho / (n0^rho + n1^rho)`. Larger `rho`
#' forces balance harder; `rho = 0` gives simple randomization and
#' `rho = 1` reproduces Wei's urn in the limit.
#'
#' @param n Number of subjects.
#' @param rho Non-negative exponent controlling the force toward
#'   balance. Default 1.
#' @return An integer vector of length `n`.
#' @references Smith RL (1984). Sequential treatment allocation
#'   using biased coin designs. *JRSS-B* 46(3):519-543.
#' @examples
#' set.seed(1)
#' alloc_smith(20, rho = 2)
#' @export
alloc_smith <- function(n, rho = 1) {
  n <- check_n(n)
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) ||
        rho < 0) {
    stop("`rho` must be a single non-negative number.",
         call. = FALSE)
  }
  out <- integer(n)
  n1 <- 0L
  n0 <- 0L
  for (i in seq_len(n)) {
    # Both counts are zero for the first subject, where the rule
    # is undefined; a fair coin is the natural completion.
    prob1 <- if (n0 == 0L && n1 == 0L) 0.5 else
      n0^rho / (n0^rho + n1^rho)
    a <- as.integer(stats::runif(1) < prob1)
    out[i] <- a
    if (a == 1L) n1 <- n1 + 1L else n0 <- n0 + 1L
  }
  out
}
