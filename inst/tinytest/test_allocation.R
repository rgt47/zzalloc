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
nn <- 4L
bb <- 2L
grid <- expand.grid(rep(list(c(0L, 1L)), nn))
adm <- grid[apply(grid, 1, function(r)
  max(abs(cumsum(2 * r - 1))) <= bb), , drop = FALSE]
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
