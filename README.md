# zzalloc

Treatment allocation schemes for clinical trials.

A unified implementation of treatment allocation procedures for
randomized clinical trials, spanning unrestricted randomization,
restricted procedures that bound imbalance in the total sample size,
stratified designs, the covariate-adaptive minimization family, and
covariate-adjusted response-adaptive designs.

Every procedure shares one sequential interface and returns
assignments in arrival order, so that procedures can be substituted
for one another in simulation studies without a `switch` statement at
each call site.

## Installation

```r
# install.packages("pak")
pak::pak("rgt47/zzalloc")
```

## Usage

Each procedure has its own function, and `allocate()` dispatches to
any of them by name:

```r
library(zzalloc)

set.seed(1)
allocate("simple", n = 10)
#>  [1] 1 1 0 0 1 0 0 0 0 1

cv <- data.frame(sex = sample(c("F", "M"), 20, TRUE))
allocate("pocock_simon", covariates = cv, p = 0.9)
#>  [1] 1 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 1 1 0 0
```

`allocation_schemes()` lists what is available:

```r
allocation_schemes()
#>  [1] "simple"            "random_allocation" "permuted_block"
#>  [4] "efron"             "wei_urn"           "smith"
#>  [7] "big_stick"         "maximal"           "stratified_block"
#> [10] "pocock_simon"      "hu_hu"             "msb"
#> [13] "cara"
```

## Available procedures

| Family | Functions |
|:-------|:----------|
| Unrestricted | `alloc_simple()`, `alloc_random_allocation()`, `alloc_permuted_block()` |
| Biased coin | `alloc_efron()`, `alloc_wei_urn()`, `alloc_smith()` |
| Restricted | `alloc_big_stick()`, `alloc_maximal()` |
| Stratified | `alloc_stratified_block()` |
| Minimization | `alloc_pocock_simon()`, `alloc_hu_hu()`, `alloc_msb()` |
| Response-adaptive | `alloc_cara()` |

## Imbalance diagnostics

`allocation_imbalance()` reports overall, marginal, and within-stratum
imbalance for an allocation. `weighted_imbalance()` applies the
response-weighted measure appropriate when covariates differ in
prognostic strength.

```r
allocation_imbalance(allocate("efron", n = 50))
```

## License

GPL-3
