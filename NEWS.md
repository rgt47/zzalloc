# zzalloc 0.1.0

First release. Seventeen exported functions providing a common
sequential interface to treatment allocation procedures for randomized
clinical trials.

## Allocation procedures

* Unrestricted randomization: `alloc_simple()`,
  `alloc_random_allocation()`, `alloc_permuted_block()`.
* Biased coin designs: `alloc_efron()`, `alloc_wei_urn()`,
  `alloc_smith()`.
* Restricted procedures bounding imbalance in the total sample size:
  `alloc_big_stick()`, `alloc_maximal()`.
* Stratified designs: `alloc_stratified_block()`.
* Covariate-adaptive minimization: `alloc_pocock_simon()`,
  `alloc_hu_hu()`, `alloc_msb()`.
* Covariate-adjusted response-adaptive designs: `alloc_cara()`.

## Interface and diagnostics

* `allocate()` provides the shared sequential entry point and
  `allocation_schemes()` enumerates the available procedures, so that
  one procedure can be substituted for another in simulation studies.
* `allocation_imbalance()` and `weighted_imbalance()` report imbalance,
  the latter using the response-weighted measure appropriate when
  covariates differ in prognostic strength.
