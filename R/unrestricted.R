#' Simple randomization
#'
#' Each subject is assigned independently by a coin flip. The
#' procedure is unpredictable and unbiased but places no bound on
#' imbalance: the difference in arm sizes has standard deviation
#' `sqrt(n)`, so imbalance grows without limit as the trial grows,
#' even though the *proportion* allocated to each arm converges.
#'
#' @param n Number of subjects.
#' @param prob Probability of assignment to treatment. Default 0.5.
#' @return An integer vector of length `n`, `0` for control and
#'   `1` for treatment, in arrival order.
#' @examples
#' set.seed(1)
#' alloc_simple(10)
#' @export
alloc_simple <- function(n, prob = 0.5) {
  n <- check_n(n)
  prob <- check_p(prob, "prob")
  as.integer(stats::runif(n) < prob)
}

#' Random allocation to a fixed arm size
#'
#' A permutation of exactly `n / 2` assignments to each arm, so
#' final imbalance is zero by construction. Unlike simple
#' randomization the assignments are not independent, and the last
#' subjects in the sequence are highly predictable once the earlier
#' ones are known.
#'
#' @param n Number of subjects. Odd `n` gives arms differing by one.
#' @return An integer vector of length `n`.
#' @examples
#' set.seed(1)
#' table(alloc_random_allocation(20))
#' @export
alloc_random_allocation <- function(n) {
  n <- check_n(n)
  n1 <- n %/% 2L
  sample(c(rep(1L, n1), rep(0L, n - n1)))
}

#' Permuted block randomization
#'
#' Within each consecutive block of `block_size` subjects the
#' assignments are a random permutation of equal numbers per arm,
#' so imbalance never exceeds half a block. Smaller blocks buy
#' tighter balance at the cost of greater predictability: with a
#' known fixed block size the final assignment in each block is
#' determined.
#'
#' @param n Number of subjects.
#' @param block_size Even block length, or a vector of permissible
#'   block lengths sampled at random for each block. Random block
#'   sizes reduce predictability.
#' @return An integer vector of length `n`.
#' @examples
#' set.seed(1)
#' alloc_permuted_block(12, block_size = 4)
#' alloc_permuted_block(12, block_size = c(2, 4, 6))
#' @export
alloc_permuted_block <- function(n, block_size = 4L) {
  n <- check_n(n)
  if (!is.numeric(block_size) || !length(block_size) ||
        anyNA(block_size)) {
    stop("`block_size` must be a numeric vector without missing ",
         "values.", call. = FALSE)
  }
  if (any(block_size < 2) || any(block_size %% 2 != 0)) {
    stop("`block_size` entries must be even and at least 2.",
         call. = FALSE)
  }
  out <- integer(0)
  while (length(out) < n) {
    b <- if (length(block_size) == 1L) block_size else
      sample(block_size, 1L)
    out <- c(out, sample(rep(c(0L, 1L), each = b / 2)))
  }
  out[seq_len(n)]
}
