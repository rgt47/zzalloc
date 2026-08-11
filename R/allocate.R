#' Allocate by a named scheme
#'
#' A single entry point to every procedure in the package, so that
#' a simulation can iterate over schemes without a switch
#' statement at each call site. Arguments specific to a scheme are
#' passed through `...`.
#'
#' @param scheme One of the values returned by
#'   [allocation_schemes()].
#' @param n Number of subjects. Required for schemes that do not
#'   take covariates; ignored, with a warning if inconsistent, for
#'   those that do.
#' @param covariates Data frame of balancing factors, required for
#'   the covariate-adaptive and stratified schemes.
#' @param ... Passed to the underlying function.
#' @return An integer vector of assignments. For `"cara"` a list,
#'   as documented in [alloc_cara()].
#' @examples
#' set.seed(1)
#' allocate("simple", n = 10)
#' cv <- data.frame(sex = sample(c("F", "M"), 20, TRUE))
#' allocate("pocock_simon", covariates = cv, p = 0.9)
#' @export
allocate <- function(scheme, n = NULL, covariates = NULL, ...) {
  scheme <- match.arg(scheme, allocation_schemes())
  needs_cov <- scheme %in% c("stratified_block", "pocock_simon",
                             "hu_hu", "msb", "cara")
  if (needs_cov) {
    if (is.null(covariates)) {
      stop("scheme '", scheme, "' requires `covariates`.",
           call. = FALSE)
    }
    if (!is.null(n) && n != nrow(covariates)) {
      stop("`n` is ", n, " but `covariates` has ",
           nrow(covariates), " rows; supply one or the other.",
           call. = FALSE)
    }
  } else {
    if (is.null(n)) {
      if (is.null(covariates)) {
        stop("scheme '", scheme, "' requires `n`.", call. = FALSE)
      }
      n <- nrow(covariates)
    }
  }
  fn <- switch(scheme,
    simple = function(...) alloc_simple(n, ...),
    random_allocation = function(...) alloc_random_allocation(n),
    permuted_block = function(...) alloc_permuted_block(n, ...),
    efron = function(...) alloc_efron(n, ...),
    wei_urn = function(...) alloc_wei_urn(n, ...),
    smith = function(...) alloc_smith(n, ...),
    big_stick = function(...) alloc_big_stick(n, ...),
    maximal = function(...) alloc_maximal(n, ...),
    stratified_block = function(...)
      alloc_stratified_block(covariates, ...),
    pocock_simon = function(...)
      alloc_pocock_simon(covariates, ...),
    hu_hu = function(...) alloc_hu_hu(covariates, ...),
    msb = function(...) alloc_msb(covariates, ...),
    cara = function(...) alloc_cara(covariates, ...)
  )
  fn(...)
}

#' Available allocation schemes
#'
#' @return A character vector of scheme names accepted by
#'   [allocate()].
#' @examples
#' allocation_schemes()
#' @export
allocation_schemes <- function() {
  c("simple", "random_allocation", "permuted_block",
    "efron", "wei_urn", "smith",
    "big_stick", "maximal",
    "stratified_block", "pocock_simon", "hu_hu", "msb",
    "cara")
}
