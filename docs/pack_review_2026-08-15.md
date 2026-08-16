# zzalloc Package Review: Gap Analysis White Paper
*2026-08-15 16:40 PDT*

Referee-grade review of the R package `zzalloc` (version 0.0.0.9000,
single commit `520e210`) at `~/prj/sfw/18-zzalloc/zzalloc`. No earlier
dated review exists in `docs/`; this is the first. Every claim below is
tagged: verified (ran it, observed the output), inspected (read the
source), inferred, or unverified.

Environment for all runs: macOS ARM64, R 4.6.1 (Homebrew), with
`R_PROFILE_USER=/dev/null` to bypass the zzcollab `.Rprofile`. The
package's `renv` machinery was therefore not active during the review.

## 1. Verdict

Not ready for CRAN submission, although it is close on mechanics. The
check itself is clean:

- `R CMD check --as-cran --no-manual`: 0 ERRORS, 0 WARNINGS, 1 NOTE
  (verified). The NOTE is 'CRAN incoming feasibility': new submission
  plus 'Version contains large components (0.0.0.9000)'.

What blocks release is not the check. Three functional defects were
verified on documented paths (a silently disabled CARA adaptation under
reserved covariate names, a numeric-overflow crash in `alloc_maximal`
at realistic trial sizes, and silent argument swallowing in
`allocate('random_allocation', ...)`), and the package has no front
door at all: no `README.md`, no `NEWS.md`, no package-level help topic,
and an empty `vignettes/` directory. The DESCRIPTION lacks
`URL`/`BugReports`, the version is a devel version, and the LICENSE
file is an MIT-style stub inconsistent with the declared GPL-3.

## 2. `R CMD check` and tooling results

### 2.1 rcmdcheck (verified)

Command:

```r
rcmdcheck::rcmdcheck('/Users/zenn/prj/sfw/18-zzalloc/zzalloc',
  args = c('--as-cran', '--no-manual'), error_on = 'never')
```

Result: 0 ERRORS, 0 WARNINGS, 1 NOTE, verbatim:

```
checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ronald G. Thomas <rgthomas@ucsd.edu>'
New submission
Version contains large components (0.0.0.9000)
```

Diagnosis and fix: set a release version (e.g. 0.1.0). The 'new
submission' half of the NOTE is unavoidable and harmless. The test
suite ran under check (`Running 'tinytest.R' [28s/89s]`, verified) and
all examples ran (`checking examples ... OK`, verified); no
`\dontrun{}` or `\donttest{}` exists anywhere (verified by grep), so
nothing was skipped. The check also noted and removed the empty
`vignettes/` directory during build (verified in the check log); the
empty directory should be removed or populated.

Suggests (`tinytest`) was installed, so nothing was silently skipped
for missing Suggests (verified).

### 2.2 DESCRIPTION and tarball hygiene (inspected)

- Version `0.0.0.9000`: not a release version. Blocker for submission.
- No `URL` or `BugReports` fields. CRAN does not require them, but
  their absence removes the standard path for bug reports.
- License: `GPL-3` in DESCRIPTION, but the root `LICENSE` file is the
  two-line `YEAR:`/`COPYRIGHT HOLDER:` MIT template stub. That stub is
  meaningless under GPL-3. Because `LICENSE` is in `.Rbuildignore` the
  tarball is internally consistent and the check passes, but the
  repository is not: delete the stub, or switch to `GPL-3 + file
  LICENSE` with a real file.
- `.Rbuildignore` correctly excludes the zzcollab scaffolding
  (`analysis/`, `docs/`, `Makefile`, `.devcontainer`, `.github`,
  `CITATION.cff`, renv artifacts). Inspected; sound.
- `CITATION.cff` exists at the root but there is no `inst/CITATION`;
  users calling `citation('zzalloc')` will get only the auto-generated
  entry. Minor.
- Imports contains only `stats`, and every cross-namespace call in
  `R/` uses `stats::` (inspected across all nine files). Sound.
- `NEWS.md`: absent. Needed for release.
- Non-ASCII characters in `R/` and DESCRIPTION: none (verified by
  grep). Trailing newlines: all `.R` files end with one (verified).
- No `.onLoad`/`.onAttach`/`zzz.R` at all (verified by directory
  listing), so no startup-message or options-restoration issues.

### 2.3 Automated tooling

- `covr::package_coverage()` (verified): 83.89% total. By file:
  imbalance.R 59.02%, utils.R 78.95%, cara.R 84.62%, biased_coin.R
  88.10%, restricted.R 88.24%, stratified.R 88.24%, allocate.R 89.47%,
  unrestricted.R 90.48%, minimization.R 92.59%. The imbalance.R gap is
  mostly unexercised error branches and the factor-covariate path of
  `weighted_imbalance()` (inferred from reading the file against the
  test suite).
- `lintr::lint_package()` (verified): 108 lints. Breakdown: 53
  object_usage (nearly all false positives; the test-file ones arise
  because lintr could not find an installed zzalloc, and the `R/`
  ones flag `n`/`covariates` captured by the closures in
  `allocate()`), 46 indentation (all in `inst/tinytest/`, which is
  indented to a different convention than lintr's default), 8 brace
  (multi-line anonymous functions without braces, mostly the
  `allocate()` dispatch table), 1 infix spacing. Nothing is a
  correctness finding; the volume is style noise worth either fixing
  or configuring away with a `.lintr` file.
- `spelling::spell_check_package()` (verified): 29 flagged words. Two
  are genuine British spellings that violate the project's US-English
  standard: 'colour' (four occurrences, `man/alloc_wei_urn.Rd`, from
  `R/biased_coin.R`; the source comment 'favouring' at
  `R/biased_coin.R:60` is a third, found by inspection) and
  'enrolment' (`man/alloc_cara.Rd:58`, from `R/cara.R`). The rest are
  proper nouns and journal names that belong in `inst/WORDLIST`.
- `urlchecker::url_check()` (verified): all URLs correct (there are
  effectively none to check).
- `codetools::checkUsagePackage('zzalloc', all = TRUE)` (verified):
  only 'parameter changed by assignment' notes (deliberate
  validate-and-coerce idiom) and one benign note on the `fn` dispatch
  variable in `allocate()`. No 'no visible binding' defects.
- `goodpractice::gp()`: NOT RUN; the package is not installed on this
  machine and installing its dependency tree was out of scope per the
  review protocol.

## 3. Functional bugs (these matter most)

### 3.1 `alloc_cara()` silently disables adaptation when a covariate is named `trt` or `resp` (verified)

`alloc_cara()` builds its working-model frame as
`cbind(covariates[...], trt = ..., resp = ...)` and the formula as
`resp ~ trt + <covariate names>` (`R/cara.R:94-99`). If the user's
covariate data frame contains a column named `trt` (or `resp`), the
bound frame has two columns with the same name; model.frame resolves
to the first, which is the covariate, so the realized assignment
never enters the model. Adaptation is silently destroyed and the
procedure quietly behaves like a clipped fair coin. Demonstration
(seed 99, n = 200, burn-in 30, strongly superior treatment,
`plogis(-1 + 3 * trt)`):

```
tilt with covariate named 'x':    0.676
tilt with covariate named 'trt':  0.535
tilt with covariate named 'resp': 0.535
```

This is the worst class of defect: a plausible input, no error, no
warning, and a wrong randomization sequence in a clinical-trial
package. Fix: reject (or rename) covariates named `trt`/`resp` at
entry, and build the formula with backticked names.

Related, same root cause: a covariate with a non-syntactic name (e.g.
`'my x'`) fails with the raw parse error `unexpected symbol` from
`as.formula()` rather than a useful message (verified).

### 3.2 `alloc_maximal()` crashes from numeric overflow at realistic n (verified)

The backward recursion stores path counts in doubles
(`R/restricted.R:96-107`). Counts grow geometrically, overflow to
`Inf`, the transition probability becomes `Inf/(Inf + Inf) = NaN`,
and the function dies with the raw message `missing value where
TRUE/FALSE needed`. Verified brackets: n = 1100 with boundary 3
works, n = 1200 fails; n = 900 with boundary 10 works, n = 1100
fails. A 1,200-subject trial is not exotic, and the failure message
gives the user no clue. Fix: normalize each row of `counts` (only the
ratio `up/(up + down)` is ever used, so per-row rescaling is exact),
or work in log space; at minimum, detect non-finite counts and error
with an informative message.

### 3.3 `allocate('random_allocation', ...)` silently swallows every extra argument (verified)

The dispatch entry is `function(...) alloc_random_allocation(n)`
(`R/allocate.R:50`), so the `...` is dropped.
`allocate('random_allocation', n = 10, block_size = 4)` returns
normally, and so would a misspelled argument intended for another
scheme. Every other dispatch entry forwards `...` and a typo is
caught downstream as an 'unused argument' error (verified with `prb =`
and `blocksize =`). Fix: forward `...` so the downstream 0-argument
signature produces the error, e.g. `function(...)
alloc_random_allocation(n, ...)`.

### 3.4 Documented behavior mismatch in `allocate()` (verified)

The `n` docs promise `n` is 'ignored, with a warning if inconsistent'
for covariate schemes (`R/allocate.R:11-13`), but the code stops with
an error (`R/allocate.R:33-37`; verified: `allocate('pocock_simon',
n = 5, covariates = cv20)` errors). The error is arguably the better
behavior; the documentation must be corrected to match it.

### 3.5 Lesser verified defects on documented paths

- `alloc_stratified_block()` validates `block_size` lazily inside the
  allocation loop (`R/stratified.R:36-44`) rather than at entry as
  `alloc_permuted_block()` does. Consequences, both verified:
  `block_size = 'a'` dies with the raw `non-numeric argument to
  binary operator`; and with a mixed vector like `c(2, 3)` the error
  arrives only when the invalid size happens to be drawn, so the
  failure point depends on the RNG state. Validate at entry with the
  same check `alloc_permuted_block()` uses.
- `allocate('pocock_simon', n = c(20, 20), covariates = cv)` fails
  with the raw coercion error `'length = 2' in coercion to
  'logical(1)'` because `n` is compared to `nrow(covariates)` before
  any validation (`R/allocate.R:33`; verified).
- `alloc_cara()` does not validate `burn_in`: `burn_in = NA` produces
  the raw `missing value where TRUE/FALSE needed`, and negative
  values are accepted silently (both verified). `alloc_msb()`
  validates its own `burn_in` correctly (verified), so the two are
  inconsistent.
- `alloc_cara()` silently accepts a `response_fn` that returns
  non-binary values (verified with a constant 0.5): the binomial GLM
  warnings are suppressed by design (`R/cara.R:103`), so a
  misspecified response function degrades adaptation without any
  signal. The docs say the response is binary; the code should
  enforce it.

## 4. Help system

Counts (verified against `NAMESPACE` and `man/`):

- 17 exports; 17 Rd topics; 17/17 have `\examples`; 17/17 have
  `\value`; 0/17 have `\seealso` or `\family`.
- Doc-vs-code mismatches found: 1 substantive (Section 3.4). For a
  sample read closely against the implementations (`alloc_simple`,
  `alloc_efron`, `alloc_maximal`, `alloc_pocock_simon`, `alloc_hu_hu`,
  `alloc_msb`, `alloc_cara`, `allocation_imbalance`,
  `weighted_imbalance`), the documented defaults, argument sets, and
  return structures match the code (inspected). Return values are
  described concretely (element names for the `alloc_cara` list, all
  four components of `allocation_imbalance`), not as bare class names.

Quality of what exists is high: every topic carries genuinely
explanatory prose about the statistical trade-offs, primary
references with full citations, and runnable seeded examples. All
examples execute under check (verified).

Gaps, in order of cost:

- No package-level topic. There is no `zzalloc-package.Rd` and no
  `@keywords internal` `_PACKAGE` roxygen block, so `?zzalloc` fails.
  This is the single cheapest fix with the highest payoff.
- No root `README.md`. Combined with the missing package topic, a new
  user has no entry point at all; discovery of
  `allocation_schemes()`/`allocate()` requires reading the index.
- No vignettes; the `vignettes/` directory is empty (verified). The
  DESCRIPTION's headline claim, that procedures 'can be substituted
  for one another in simulation studies', is exactly the workflow a
  vignette should demonstrate end to end (compare several schemes on
  `allocation_imbalance()` and `weighted_imbalance()` across
  replicates). No export is demonstrated in any vignette because
  there are none.
- No cross-references. The prose does link related topics inline
  (e.g. `alloc_hu_hu` to `alloc_pocock_simon`, `allocate` to
  `allocation_schemes`; inspected), but with no `@family alloc`
  block a user landing on `alloc_efron` gets no systematic route to
  the other twelve schemes. One `@family` tag across the `alloc_*`
  set would fix this.
- `NEWS.md` absent.

## 5. User interface

### 5.1 First-use walkthrough (verified)

With no README and no vignette, the walkthrough must start at the help
index. From `?allocate`'s example, verbatim:

```r
set.seed(1)
allocate('simple', n = 10)
#  [1] 1 1 0 0 1 0 0 0 0 1
cv <- data.frame(sex = sample(c('F', 'M'), 20, TRUE))
allocate('pocock_simon', covariates = cv, p = 0.9)
#  [1] 1 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 1 1 0 0
```

Two lines to a first real result once you have found `?allocate`;
both ran exactly as documented. The friction is entirely upstream of
the first call: nothing tells you `allocate()` is the front door. The
finding is the missing front door (Section 4), not the API, which is
compact and learnable.

### 5.2 Findings

Strengths, tersely: consistent `alloc_` prefix (good tab completion);
data-first argument order throughout (`n` or `covariates` first;
pipe-friendly); uniform 0/1 integer coding in arrival order across
all schemes; error messages generally name the argument, state the
constraint, and use `call. = FALSE`; the export surface is minimal
(17 exports, no leaked internals; helpers are `@noRd`).

Weaknesses:

- Silent surprises: Sections 3.1, 3.3, and 3.5 list four verified
  ones. Additionally, `allocate('simple', covariates = cv)` silently
  takes `n` from `nrow(covariates)` and ignores the covariate values
  (verified); this fallback is undocumented, and accepting covariates
  for a scheme that cannot use them deserves at least a warning.
- Argument-name inconsistency: `alloc_simple()` calls its probability
  `prob` (P(treatment)) while the six schemes in the biased-coin and
  minimization families call theirs `p` (P(deficient arm)). The
  concepts differ, so distinct names are defensible, but the near
  collision is a trap: `alloc_simple(10, p = 0.9)` partial-matches
  `prob` silently (verified), so a user who habitually types `p =`
  gets a different meaning with no error. Renaming `prob` to
  something non-prefix-colliding (e.g. `p_treat`), or vice versa, is
  cheapest now, before the first release.
- Type stability of `allocate()`: 12 schemes return an integer
  vector, `'cara'` returns a list (both verified). It is documented,
  but it forces a branch in exactly the scheme-agnostic simulation
  loop the function exists to serve. Consider returning a uniform
  structure (or having `alloc_cara()` return the assignment vector
  with the responses as an attribute or separate accessor). API
  changes are cheapest before the first release.
- Print/format methods: none exist; all returns are bare vectors and
  lists. Acceptable at this scale, and `allocation_imbalance()`'s
  list prints legibly (verified), but a 5,000-subject assignment
  vector floods the console. Low priority.
- Defaults are sound and documented with their provenance (Efron's
  2/3, Pocock-Simon's `sd` measure, RSIHR target for CARA;
  inspected). `alloc_msb()`'s `threshold = 0.3` is unusual but the
  docs justify it explicitly. No cross-function default
  contradictions found.

## 6. Coding practices

Checked and found sound (all inspected, across all nine `R/` files
unless noted):

- Assignment exclusively `<-`; no `T`/`F`; `seq_len`/`seq_along`
  everywhere (no `1:n`); `vapply` rather than `sapply`; `drop =
  FALSE` on every data-frame row subset in `cara.R`; no
  `library()`/`require()` in `R/`; no `set.seed()` in package code
  (only in examples and tests; verified by grep); no
  `options()`/`par()` calls at all; no S3 classes, so no
  registration or dispatch issues arise.
- Validation helpers (`check_n`, `check_p`, `check_covariates`,
  `resolve_weights`, `check_trt`) are centralized and produce
  named-argument messages. `check_covariates()` rejecting continuous
  columns with an explanation of why (`R/utils.R:56-61`) is a model
  error message.
- RNG: all draws go through `stats::runif`/`sample`; the package
  makes no RNG-discipline claims and sets no seeds. Sound.
- Tests are behavioral, not structural: they verify imbalance bounds
  for `big_stick`/`maximal`/`permuted_block`, exact uniformity of
  the maximal procedure against an independently enumerated
  admissible set (a genuinely strong test), the balance ordering of
  minimization versus simple randomization, the weighted-minimization
  trade-off in both directions, and CARA's adaptive tilt
  (`inst/tinytest/test_allocation.R`; inspected, and verified to pass
  under `R CMD check`, not merely `load_all()`). No `:::` calls in
  tests (verified by grep), so nothing passes only under `load_all()`.

Findings:

- `R/restricted.R:96-107`: unguarded double overflow in the maximal
  procedure's count recursion (Section 3.2).
- `R/cara.R:94-99`: formula built by pasting unescaped covariate
  names into a frame that also injects `trt`/`resp` (Section 3.1).
- `R/allocate.R:50`: dropped `...` (Section 3.3).
- `R/stratified.R:36-44`: lazy `block_size` validation (Section 3.5).
- `R/allocate.R:33`: unvalidated `n` compared before coercion
  (Section 3.5).
- `R/unrestricted.R:70-77` (`alloc_permuted_block`): grows `out` by
  `c()` inside a while loop. Harmless at trial scale (at most
  `n/2 + max block` iterations) but the one performance smell
  present; preallocation is trivial. Inspected, not benchmarked.
- Untested code paths (from the covr report, verified): most error
  branches in `imbalance.R` and `utils.R`, the factor-covariate path
  of `weighted_imbalance()`, `alloc_cara()`'s degenerate-fit
  fallbacks, and the `measure = 'diff'`/`'squared'` variants, which
  no test exercises at all.
- Test-suite runtime: 28s CPU/89s elapsed under check on this
  machine, dominated by the 12,000-replicate uniformity test.
  Acceptable but near the level where CRAN may complain on slow
  platforms; consider trimming replicates or gating the heavy test
  behind `at_home()`.

## 7. Prioritized checklist

(a) CRAN blockers

- Set a release version in DESCRIPTION (0.0.0.9000 is flagged by the
  incoming check).
- Resolve the LICENSE inconsistency: delete the MIT-template stub or
  adopt `GPL-3 + file LICENSE` with a real file.
- Add `NEWS.md` and a root `README.md`.
- Remove or populate the empty `vignettes/` directory.
- Add `URL` and `BugReports` to DESCRIPTION.
- Run win-builder and R-hub (R-devel) before submitting; neither was
  run here.

(b) Bugs to fix before anyone depends on the behavior

- `alloc_cara()`: reject or safely handle covariates named
  `trt`/`resp`, and handle non-syntactic names (Section 3.1).
- `alloc_maximal()`: rescale the count recursion so n in the
  thousands works, or fail informatively (Section 3.2).
- `allocate('random_allocation', ...)`: forward `...` so typos error
  (Section 3.3).
- `alloc_stratified_block()`: validate `block_size` at entry
  (Section 3.5).
- `allocate()`: validate `n` before comparing to `nrow(covariates)`;
  fix the `n` docs ('warning' vs actual error, Section 3.4); warn or
  error when covariates are supplied to a scheme that ignores them.
- `alloc_cara()`: validate `burn_in`; enforce binary responses.

(c) Documentation completion

- Add a `_PACKAGE` roxygen block so `?zzalloc` works.
- Add `@family` tags across the `alloc_*` schemes and the imbalance
  pair.
- Write a getting-started/simulation vignette demonstrating the
  scheme-substitution workflow plus the imbalance diagnostics.
- Fix US-English drift: 'colour' (x4) and 'favouring' in
  `R/biased_coin.R`, 'enrolment' in `R/cara.R`; add the proper nouns
  from the spell check to `inst/WORDLIST`.
- Move citation metadata into `inst/CITATION` (keep `CITATION.cff`
  for GitHub if desired).

(d) Design decisions best made before first release

- The `prob` (alloc_simple) versus `p` (everything else) near
  collision, given verified silent partial matching (Section 5.2).
- Whether `allocate()` should be type-stable across all schemes,
  including `'cara'` (Section 5.2).
- Whether assignments deserve a lightweight class with a `print`
  method, which would also give `allocation_imbalance()` a natural
  `summary` home. Optional.
- Test coverage for the untested `measure` variants and error
  branches (imbalance.R at 59%).

## 8. Not evaluated

- Platforms: only macOS ARM64 / R 4.6.1 was used. NOT checked:
  win-builder, R-hub, R-devel, any Linux, any older R (the DESCRIPTION
  claims R >= 4.1.0; that floor is unverified).
- `goodpractice::gp()`: not installed, not run.
- The zzcollab container workflow (`Dockerfile`, `make r`,
  `renv`): not exercised; the review deliberately bypassed the
  project `.Rprofile`.
- The `--no-manual` flag means PDF manual construction (and any LaTeX
  issues in it) was not checked.
- Statistical fidelity of each procedure to its cited reference was
  verified only where the test suite does so (maximal-procedure
  uniformity, imbalance bounds, degenerate-p behavior); the
  Pocock-Simon, Hu-Hu, MSB, Wei, and Smith implementations were read
  against their documented descriptions (inspected) but not validated
  against published numerical results.
- Long-run distributional properties of `alloc_cara()` beyond the
  single verified tilt comparison.
- Performance beyond reading for smells; nothing was benchmarked.
- Concurrent-review interference note: this machine was running R CMD
  checks for several sibling packages during the review. The zzalloc
  check was rerun with an isolated log and check directory, and the
  counts reported in Section 1 come from that isolated run's parsed
  rcmdcheck object, not from any shared log.
