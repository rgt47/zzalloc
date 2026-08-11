#' Stratified permuted block randomization
#'
#' A separate permuted block sequence is maintained within each
#' cell of the cross-classification of the balancing factors, so
#' balance is enforced jointly rather than marginally. This is the
#' strongest form of balance available, and its weakness is
#' arithmetic: the number of cells is the product of the numbers of
#' levels, so with more than a few factors most cells hold a
#' handful of subjects and each contributes its own partial block,
#' at which point marginal balance can be worse than under
#' procedures that target the margins directly.
#'
#' @param covariates A data frame of discrete balancing factors,
#'   one row per subject in arrival order.
#' @param block_size Even block length, or a vector of permissible
#'   lengths sampled per block within each stratum.
#' @return An integer vector with one assignment per row of
#'   `covariates`, in arrival order.
#' @examples
#' set.seed(1)
#' cv <- data.frame(sex = rep(c("F", "M"), each = 10),
#'                  stage = rep(c("I", "II"), 10))
#' alloc_stratified_block(cv, block_size = 4)
#' @export
alloc_stratified_block <- function(covariates, block_size = 4L) {
  covariates <- check_covariates(covariates)
  n <- nrow(covariates)
  stratum <- interaction(covariates, drop = FALSE, sep = "|")
  out <- integer(n)
  # One independent block sequence per stratum, consumed in
  # arrival order. Queues are generated lazily so that strata
  # never visited cost nothing.
  queues <- list()
  for (i in seq_len(n)) {
    s <- as.character(stratum[i])
    if (is.null(queues[[s]]) || !length(queues[[s]])) {
      b <- if (length(block_size) == 1L) block_size else
        sample(block_size, 1L)
      if (b < 2 || b %% 2 != 0) {
        stop("`block_size` entries must be even and at least 2.",
             call. = FALSE)
      }
      queues[[s]] <- sample(rep(c(0L, 1L), each = b / 2))
    }
    out[i] <- queues[[s]][1L]
    queues[[s]] <- queues[[s]][-1L]
  }
  out
}
