#' Big stick design
#'
#' Assignments are a fair coin flip until the imbalance reaches
#' `boundary`, at which point the next subject is forced to the
#' deficient arm. The procedure is unpredictable everywhere except
#' at the boundary, which is where it differs from the biased-coin
#' family: those tilt the coin at every step, while the big stick
#' pays for its bound only when the bound binds.
#'
#' @param n Number of subjects.
#' @param boundary Maximum permitted absolute imbalance. Must be a
#'   positive whole number.
#' @return An integer vector of length `n`.
#' @references Soares JF, Wu CFJ (1983). Some restricted
#'   randomization rules in sequential designs. *Communications in
#'   Statistics* 12(17):2017-2034.
#' @examples
#' set.seed(1)
#' x <- alloc_big_stick(50, boundary = 3)
#' max(abs(cumsum(2 * x - 1)))
#' @export
alloc_big_stick <- function(n, boundary = 3L) {
  n <- check_n(n)
  if (length(boundary) != 1L || is.na(boundary) ||
        !is.numeric(boundary) || boundary < 1 ||
        boundary != floor(boundary)) {
    stop("`boundary` must be a single positive whole number.",
         call. = FALSE)
  }
  out <- integer(n)
  d <- 0L
  for (i in seq_len(n)) {
    a <- if (d >= boundary) {
      0L
    } else if (d <= -boundary) {
      1L
    } else {
      as.integer(stats::runif(1) < 0.5)
    }
    out[i] <- a
    d <- d + if (a == 1L) 1L else -1L
  }
  out
}

#' Maximal procedure
#'
#' Draws uniformly from the set of all assignment sequences whose
#' running imbalance never exceeds `boundary`. Among procedures
#' respecting that bound this is the most random one available, so
#' it minimizes the selection bias that arises when an
#' investigator can guess upcoming assignments. Berger and
#' colleagues proposed it for exactly that reason.
#'
#' The sequence is built one subject at a time, with the
#' probability of each assignment proportional to the number of
#' valid completions it leaves. Those counts are obtained by a
#' backward recursion over (subjects remaining, current
#' imbalance), which is what makes the draw exactly uniform rather
#' than merely bounded.
#'
#' @param n Number of subjects.
#' @param boundary Maximum permitted absolute imbalance.
#' @return An integer vector of length `n`.
#' @references Berger VW, Ivanova A, Knoll MD (2003). Minimizing
#'   predictability while retaining balance through the use of less
#'   restrictive randomization procedures. *Statistics in Medicine*
#'   22(19):3017-3028.
#' @examples
#' set.seed(1)
#' x <- alloc_maximal(20, boundary = 2)
#' max(abs(cumsum(2 * x - 1)))
#' @export
alloc_maximal <- function(n, boundary = 2L) {
  n <- check_n(n)
  if (length(boundary) != 1L || is.na(boundary) ||
        !is.numeric(boundary) || boundary < 1 ||
        boundary != floor(boundary)) {
    stop("`boundary` must be a single positive whole number.",
         call. = FALSE)
  }
  boundary <- as.integer(boundary)
  states <- seq.int(-boundary, boundary)
  idx <- function(d) d + boundary + 1L

  # counts[k, ] holds, for each reachable imbalance, the number of
  # ways to allocate the remaining k subjects without ever leaving
  # [-boundary, boundary]. Row 1 is "no subjects remaining".
  counts <- matrix(0, nrow = n + 1L, ncol = length(states))
  counts[1L, ] <- 1
  for (k in seq_len(n)) {
    for (d in states) {
      tot <- 0
      if (d + 1L <= boundary) tot <- tot + counts[k, idx(d + 1L)]
      if (d - 1L >= -boundary) tot <- tot + counts[k, idx(d - 1L)]
      counts[k + 1L, idx(d)] <- tot
    }
  }
  if (counts[n + 1L, idx(0L)] == 0) {
    stop("no sequence of length ", n, " keeps the imbalance ",
         "within ", boundary, ".", call. = FALSE)
  }

  out <- integer(n)
  d <- 0L
  for (i in seq_len(n)) {
    remaining <- n - i
    up <- if (d + 1L <= boundary)
      counts[remaining + 1L, idx(d + 1L)] else 0
    down <- if (d - 1L >= -boundary)
      counts[remaining + 1L, idx(d - 1L)] else 0
    a <- as.integer(stats::runif(1) < up / (up + down))
    out[i] <- a
    d <- d + if (a == 1L) 1L else -1L
  }
  out
}
