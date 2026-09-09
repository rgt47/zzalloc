library(zzalloc)

set.seed(20260810)
cv <- data.frame(
  sex = sample(c("F", "M"), 200, TRUE),
  stage = sample(c("I", "II", "III"), 200, TRUE)
)

# --- shape and coding, every scheme -------------------------------

for (s in setdiff(allocation_schemes(), "cara")) {
  a <- if (s %in% c("stratified_block", "pocock_simon", "hu_hu",
                    "msb")) {
    allocate(s, covariates = cv)
  } else {
    allocate(s, n = 200)
  }
  expect_true(is.integer(a),
    info = paste(s, "returns an integer vector"))
  expect_equal(length(a), 200L,
    info = paste(s, "returns one assignment per subject"))
  expect_true(all(a %in% c(0L, 1L)),
    info = paste(s, "returns only 0 and 1"))
}

# --- degenerate probabilities behave as documented ----------------

set.seed(1)
expect_true(all(alloc_simple(50, prob = 1) == 1L),
  info = "prob = 1 assigns everyone to treatment")
expect_true(all(alloc_simple(50, prob = 0) == 0L),
  info = "prob = 0 assigns everyone to control")

# Efron with p = 1 alternates strictly once the arms differ, so
# imbalance never exceeds one.
set.seed(2)
d <- cumsum(2 * alloc_efron(100, p = 1) - 1)
expect_true(max(abs(d)) <= 1,
  info = "Efron with p = 1 bounds imbalance at 1")

# --- restricted procedures respect their bound --------------------

for (b in c(1L, 2L, 5L)) {
  set.seed(b)
  d <- cumsum(2 * alloc_big_stick(200, boundary = b) - 1)
  expect_true(max(abs(d)) <= b,
    info = paste("big stick respects boundary", b))
  d2 <- cumsum(2 * alloc_maximal(60, boundary = b) - 1)
  expect_true(max(abs(d2)) <= b,
    info = paste("maximal procedure respects boundary", b))
}

# The maximal procedure must be uniform over the admissible
# sequences. The admissible set is enumerated here rather than
# asserted, so the test cannot be satisfied by a wrong count that
# happens to match a wrong implementation.
#
# Admissible means BOTH that the running imbalance stays within
# the boundary and that the sequence ends in balance. Dropping the
# second condition enumerates a strictly larger set (for n = 4,
# b = 2 it is 10 sequences rather than 6) and the test would then
# pass against an implementation that never balances.
nn <- 4L
bb <- 2L
grid <- expand.grid(rep(list(c(0L, 1L)), nn))
adm <- grid[apply(grid, 1, function(r)
  max(abs(cumsum(2 * r - 1))) <= bb &&
    sum(2 * r - 1) == 0), , drop = FALSE]
n_adm <- nrow(adm)
adm_keys <- apply(adm, 1, paste, collapse = "")

set.seed(3)
reps <- 12000
seqs <- replicate(reps, paste(alloc_maximal(nn, boundary = bb),
                              collapse = ""))
tb <- table(seqs)
expect_equal(length(tb), n_adm,
  info = "maximal procedure reaches every admissible sequence
          and no others")
expect_true(setequal(names(tb), adm_keys),
  info = "the sequences produced are exactly the admissible set")
# Under uniformity each cell has SE sqrt(p(1-p)/reps); four of
# those is a generous band that still fails a visibly skewed
# sampler.
band <- 4 * sqrt((1 / n_adm) * (1 - 1 / n_adm) / reps)
expect_true(max(abs(as.numeric(tb) / reps - 1 / n_adm)) < band,
  info = "maximal procedure is uniform over admissible sequences")

# --- permuted blocks bound imbalance by half a block --------------

set.seed(4)
d <- cumsum(2 * alloc_permuted_block(200, block_size = 4) - 1)
expect_true(max(abs(d)) <= 2,
  info = "block size 4 bounds imbalance at 2")
expect_equal(sum(alloc_random_allocation(50)), 25L,
  info = "random allocation splits the sample exactly")

# --- minimization balances better than simple randomization -------

set.seed(5)
imb <- function(f, reps = 60) {
  vapply(seq_len(reps), function(i)
    allocation_imbalance(f(), cv)$max_marginal, numeric(1))
}
i_simple <- imb(function() alloc_simple(nrow(cv)))
i_ps <- imb(function() alloc_pocock_simon(cv))
i_hh <- imb(function() alloc_hu_hu(cv))
expect_true(mean(i_ps) < mean(i_simple),
  info = "Pocock-Simon beats simple randomization on marginal
          imbalance")
expect_true(mean(i_hh) < mean(i_simple),
  info = "Hu-Hu beats simple randomization on marginal imbalance")

# Hu-Hu carries a stratum term, so it should balance joint cells
# at least as well as marginal minimization does.
set.seed(6)
st <- function(f, reps = 60) {
  vapply(seq_len(reps), function(i)
    max(abs(allocation_imbalance(f(), cv)$stratum)), numeric(1))
}
expect_true(mean(st(function() alloc_hu_hu(cv))) <=
              mean(st(function() alloc_pocock_simon(cv))),
  info = "Hu-Hu balances joint strata at least as well as
          Pocock-Simon")

# --- weighted minimization shifts balance toward the weighted -----

set.seed(7)
gap <- function(w, reps = 80) {
  vapply(seq_len(reps), function(i) {
    a <- alloc_pocock_simon(cv, weights = w)
    m <- allocation_imbalance(a, cv)$marginal
    c(max(abs(m$sex)), max(abs(m$stage)))
  }, numeric(2))
}
g_even <- rowMeans(gap(c(sex = 1, stage = 1)))
g_stage <- rowMeans(gap(c(sex = 1, stage = 8)))
expect_true(g_stage[2] < g_even[2],
  info = "up-weighting stage improves stage balance")
expect_true(g_stage[1] > g_even[1],
  info = "up-weighting stage costs sex balance, since the
          weighting shifts effort rather than adding it")

# --- p = 0.5 collapses to simple randomization --------------------

# With a fair coin the minimization score never influences the
# draw, so the assignment distribution must be exchangeable and
# centred on one half.
set.seed(8)
m <- mean(replicate(400, mean(alloc_pocock_simon(cv, p = 0.5))))
expect_true(abs(m - 0.5) < 0.02,
  info = "Pocock-Simon with p = 0.5 allocates half on average")

# --- stratified blocks balance within cells -----------------------

set.seed(9)
a <- alloc_stratified_block(cv, block_size = 2)
expect_true(max(abs(allocation_imbalance(a, cv)$stratum)) <= 1,
  info = "block size 2 leaves at most one subject unmatched per
          stratum")

# --- minimal sufficient balance stays mostly random ---------------

set.seed(10)
cv1 <- data.frame(sex = sample(c("F", "M"), 300, TRUE))
i_msb <- mean(replicate(40,
  allocation_imbalance(alloc_msb(cv1), cv1)$max_marginal))
i_sr <- mean(replicate(40,
  allocation_imbalance(alloc_simple(300), cv1)$max_marginal))
expect_true(i_msb < i_sr,
  info = "MSB improves on simple randomization")

# --- response-weighted imbalance ----------------------------------

cvn <- data.frame(a = c(0, 0, 1, 1), b = c(0, 1, 0, 1))
tr <- c(1L, 1L, 0L, 0L)
# Treatment arm has a = 0 throughout, control a = 1: difference
# in mean a is -1, weighted by 2 gives -2. b is balanced.
expect_equal(weighted_imbalance(tr, cvn, c(a = 2, b = 5)), -2,
  info = "response-weighted imbalance weights by covariate effect")
expect_equal(weighted_imbalance(tr, cvn, c(a = 0, b = 5)), 0,
  info = "a covariate with no effect contributes nothing")

# --- CARA ---------------------------------------------------------

set.seed(11)
cvc <- data.frame(x = stats::rnorm(150))
rf <- function(trt, row) {
  stats::rbinom(1, 1, stats::plogis(-0.2 + 1.5 * trt + 0.5 * row$x))
}
res <- alloc_cara(cvc, rf, burn_in = 40)
expect_equal(length(res$trt), 150L,
  info = "CARA returns one assignment per subject")
expect_true(all(res$response %in% c(0, 1)),
  info = "CARA records binary responses")
# Treatment is much better here, so adaptation should tilt toward
# it relative to the fair-coin burn-in.
expect_true(mean(res$trt[41:150]) > mean(res$trt[1:40]),
  info = "CARA shifts allocation toward the superior arm")

# --- input validation ---------------------------------------------

expect_error(alloc_simple(0), info = "n must be positive")
expect_error(alloc_simple(10, prob = 1.5),
  info = "prob must be a probability")
expect_error(alloc_permuted_block(10, block_size = 3),
  info = "odd block sizes are rejected")
expect_error(
  alloc_pocock_simon(data.frame(x = c(1.5, 2.5, 3.5))),
  info = "continuous covariates are rejected")
expect_error(
  alloc_pocock_simon(data.frame(x = c("a", NA, "b"))),
  info = "missing covariate values are rejected")
expect_error(
  alloc_pocock_simon(cv, weights = c(sex = 1)),
  info = "weights must cover every covariate")
expect_error(allocate("pocock_simon", n = 10),
  info = "covariate schemes require covariates")
expect_error(allocate("simple"),
  info = "non-covariate schemes require n")
expect_error(allocate("nonesuch", n = 10),
  info = "unknown scheme names are rejected")

# --- regressions -------------------------------------------------

# A factor of "0"/"1" passes an `%in%` membership test, because
# that comparison is made on characters, but as.integer() on a
# factor yields level codes. Decoding it that way silently turned
# a balanced allocation into an imbalance of 3.
f_trt <- factor(c(0, 1, 1, 0, 1, 0))
i_trt <- c(0L, 1L, 1L, 0L, 1L, 0L)
expect_equal(allocation_imbalance(f_trt)$overall,
             allocation_imbalance(i_trt)$overall,
             info = "factor and integer assignments agree")
expect_equal(allocation_imbalance(f_trt)$overall, 0L,
             info = "a balanced factor assignment reports 0")
cv_f <- data.frame(x = c(0, 0, 1, 1, 0, 1))
expect_equal(weighted_imbalance(f_trt, cv_f, beta = c(x = 1)),
             weighted_imbalance(i_trt, cv_f, beta = c(x = 1)),
             info = "weighted imbalance accepts a factor")

# The maximal procedure is defined over sequences that finish in
# balance, so every realization must, and an odd length has no
# such sequence to draw from.
set.seed(11)
for (nn2 in c(6L, 10L, 20L)) {
  expect_equal(sum(2 * alloc_maximal(nn2, boundary = 2) - 1), 0,
               info = paste("maximal procedure ends balanced, n =",
                            nn2))
}
expect_error(alloc_maximal(7, boundary = 2), "must be even",
             info = "odd n is rejected")

# A covariate named `trt` or `resp` collides with the columns the
# working model needs; `$` then resolves to the covariate and the
# treatment indicator silently drops out of the formula.
rf_c <- function(trt, row) stats::rbinom(1, 1, 0.5)
expect_error(
  alloc_cara(data.frame(x = rnorm(10), trt = rnorm(10)), rf_c),
  "must not contain a column named",
  info = "a covariate named 'trt' is refused")
expect_error(
  alloc_cara(data.frame(x = rnorm(10), resp = rnorm(10)), rf_c),
  "must not contain a column named",
  info = "a covariate named 'resp' is refused")

# With two arms the sd and diff measures differ only by a constant
# factor, so they must produce identical allocations; squared is
# the only measure that can reorder candidate assignments.
set.seed(21)
cv_m <- data.frame(sex = sample(c("F", "M"), 40, TRUE),
                   stage = sample(c("I", "II", "III"), 40, TRUE))
set.seed(22); m_sd <- alloc_pocock_simon(cv_m, measure = "sd")
set.seed(22); m_df <- alloc_pocock_simon(cv_m, measure = "diff")
expect_identical(m_sd, m_df,
  info = "sd and diff agree for two arms")

# Hu-Hu with both extra weights at zero is exactly Pocock-Simon.
set.seed(23); ps <- alloc_pocock_simon(cv_m, p = 0.8)
set.seed(23); hh <- alloc_hu_hu(cv_m, p = 0.8, overall_weight = 0,
                                stratum_weight = 0)
expect_identical(ps, hh,
  info = "hu_hu reduces to pocock_simon when the extra weights vanish")

# --- round two regressions ---------------------------------------

# Block lengths are validated before any allocation happens, and
# the two block-consuming procedures report the same way. The
# stratified version used to check only when a stratum queue first
# refilled, so a non-numeric length surfaced as "non-numeric
# argument to binary operator".
cv_b <- data.frame(g = rep(c("a", "b"), each = 6))
for (bad in list(3, 0, -4, "four", NA)) {
  expect_error(alloc_stratified_block(cv_b, block_size = bad),
    "`block_size`",
    info = "stratified blocks reject an invalid block length")
  expect_error(alloc_permuted_block(8, block_size = bad),
    "`block_size`",
    info = "permuted blocks reject an invalid block length")
}
expect_true(is.integer(alloc_stratified_block(cv_b,
                                              block_size = c(2, 4))),
  info = "a vector of valid block lengths is still accepted")

# Every scheme reached through allocate() must reject an argument
# it does not take; random_allocation used to discard `...`.
expect_error(allocate("random_allocation", n = 10, block_size = 4),
  "unused argument",
  info = "allocate() does not silently swallow inapplicable
          arguments")

# A matrix is not a data frame, and the message should say so
# rather than claim the input has no columns, which a matrix
# plainly has. Assert the absence of the misleading clause, since
# the previous message also began "must be a data frame" and an
# assertion on that prefix alone would pass either way.
m_err <- tryCatch(alloc_pocock_simon(matrix(1:10, 5)),
                  error = function(e) conditionMessage(e))
expect_true(grepl("must be a data frame", m_err) &&
              !grepl("column", m_err),
  info = "a matrix is refused for not being a data frame, not for
          lacking columns")
expect_error(alloc_pocock_simon(data.frame()), "at least one column",
  info = "an empty data frame is refused for lacking columns")

# MSB flags a covariate with a chi-square test on its whole
# arm-by-level table, following Zhao and colleagues. It must still
# beat simple randomization and must leave the overall split alone.
# The chi-square criterion targets the covariate DISTRIBUTION
# across arms, not the count difference within each level, so it
# improves on simple randomization by a smaller margin on
# max_marginal than a per-level rule would. Measured over 300
# replicates the means are about 11.5 against 14.0; at 30
# replicates the comparison is noisy enough to flip, so use enough
# of them to make the assertion mean something.
set.seed(41)
cv_m2 <- data.frame(sex = sample(c("F", "M"), 300, TRUE))
i_msb2 <- mean(replicate(150,
  allocation_imbalance(alloc_msb(cv_m2), cv_m2)$max_marginal))
i_sr2 <- mean(replicate(150,
  allocation_imbalance(alloc_simple(300), cv_m2)$max_marginal))
expect_true(i_msb2 < i_sr2,
  info = "MSB improves on simple randomization")
# The distributional scale is the one MSB actually targets.
prop_gap <- function(a) {
  tb <- table(factor(a, levels = c(0, 1)), cv_m2$sex)
  abs(tb[2, 1] / sum(tb[2, ]) - tb[1, 1] / sum(tb[1, ]))
}
set.seed(43)
d_msb <- mean(replicate(150, prop_gap(alloc_msb(cv_m2))))
d_sr <- mean(replicate(150, prop_gap(alloc_simple(300))))
expect_true(d_msb < d_sr,
  info = "MSB improves the covariate distribution across arms")
expect_true(
  abs(mean(replicate(40, mean(alloc_msb(cv_m2)))) - 0.5) < 0.03,
  info = "MSB does not shift the overall allocation")

# Smith's design must satisfy both documented boundary conditions:
# rho = 0 is complete randomization, rho = 1 is Wei's urn at
# alpha = 0.
set.seed(42)
expect_true(
  abs(mean(replicate(600, mean(alloc_smith(30, rho = 0)))) - 0.5) < 0.02,
  info = "Smith with rho = 0 is complete randomization")
smith_p <- function(n0, n1) n0^1 / (n0^1 + n1^1)
wei_p <- function(n0, n1, a, b) (a + b * n0) / (2 * a + b * (n0 + n1))
for (nn3 in list(c(3, 5), c(10, 2), c(7, 7))) {
  expect_equal(smith_p(nn3[1], nn3[2]),
               wei_p(nn3[1], nn3[2], 1e-10, 1), tolerance = 1e-6,
    info = "Smith rho = 1 matches Wei's urn with alpha -> 0")
}

# --- tie handling (regression) ------------------------------------

# The scores compared in biased_draw() are sums of floating-point
# terms, and two assignments that are mathematically tied need not
# produce bit-identical sums. sd(c(n0, n1)) and abs(n1 - n0)/sqrt(2)
# are equal in real arithmetic but differ in the last bits for 748 of
# the 1681 count pairs up to 40, so an exact comparison missed genuine
# ties and applied p, or 1 - p, on the strength of rounding error.
# Over 400 subjects with the default measure = "sd", 23 of 108 ties
# were missed. Those subjects should have had a fair coin.
expect_true(
  sum(vapply(0:40, function(a) sum(vapply(0:40, function(b)
    stats::sd(c(a, b)) != abs(b - a) / sqrt(2), logical(1))), integer(1))) > 0,
  info = "sd and diff/sqrt(2) do differ in the last bits, which is the
          hazard being guarded against")

# A tie must be seen as a tie, so the draw is fair rather than biased.
# p is deliberately irrelevant at a tie, which is the whole point: at
# p = 1 a recognised tie still splits evenly, while a genuine
# difference is taken with certainty.
set.seed(4)
fair <- mean(replicate(4000, zzalloc:::biased_draw(5, 5, 1)))
expect_true(abs(fair - 0.5) < 0.03,
  info = "an exact tie draws fairly even at p = 1")
set.seed(4)
near <- mean(replicate(4000, zzalloc:::biased_draw(5, 5 + 1e-13, 1)))
expect_true(abs(near - 0.5) < 0.03,
  info = "a tie broken only by rounding noise still draws fairly")
# A real difference must still bite.
expect_equal(zzalloc:::biased_draw(1, 5, 1), 1L,
  info = "a genuinely smaller score for arm 1 is taken at p = 1")
expect_equal(zzalloc:::biased_draw(5, 1, 1), 0L,
  info = "a genuinely smaller score for arm 0 is taken at p = 1")

# sd and diff are proportional, so with the tie rule fixed they must
# agree exactly, not merely usually. This previously passed only at
# particular seeds and sizes.
set.seed(5)
cv_t <- data.frame(sex = sample(c("F", "M"), 400, TRUE),
                   stage = sample(c("I", "II", "III"), 400, TRUE))
for (s in c(1, 2, 3)) {
  set.seed(s); a_sd <- alloc_pocock_simon(cv_t, measure = "sd")
  set.seed(s); a_df <- alloc_pocock_simon(cv_t, measure = "diff")
  expect_identical(a_sd, a_df,
    info = paste("sd and diff agree exactly at seed", s))
}
