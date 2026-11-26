# SWIQ

**S**ample-**W**ise **I**mportance **Q**uantification (SWIQ) is a random-forest-based test for association between markers (e.g. SNP genotypes) and a quantitative trait. This repository contains the R implementation used in the accompanying paper.

## How it works

SWIQ fits a regression forest with [ranger](https://github.com/imbs-hl/ranger) and, for each marker, compares each out-of-bag sample's prediction error under its true allele with the error blended across all possible alleles, using only the trees that split on that marker. A marker's per-sample importances are then tested with a one-sided binomial (sign) test, giving one p-value per marker.

Optional population-structure covariates are forced into every tree (`always.split.variables`) but are not themselves tested.

## Files

- `scripts/swiq.R`: pure R implementation.
- `scripts/swiq_rcpp.R`: same interface, with the per-marker computation done in C++.
- `scripts/swiq_importances.cpp`: C++/OpenMP code used by `swiq_rcpp.R`. It reproduces the R version's arithmetic exactly, so both give identical p-values.

## Requirements

R packages `ranger`, `abind` and `parallel`, plus `Rcpp` (and a C++ compiler) for the faster version.

## Usage

```r
# Pure R version
source("scripts/swiq.R")

# Or the faster Rcpp version (compile once per session)
Rcpp::sourceCpp("scripts/swiq_importances.cpp")
source("scripts/swiq_rcpp.R")

p <- swiq(y = phenotype,   # numeric trait vector
          x = genotypes,   # samples x markers (data.frame or matrix)
          num.trees = 10000,
          num.threads = 4,
          pop.str = pcs)   # optional, e.g. principal components

# p: one p-value per column of x
```

Both scripts define a function called `swiq()`, so source only one of them.

### Main arguments

| Argument | Default | Description |
|---|---|---|
| `num.trees` | 10000 | Trees per forest |
| `num.threads` | 1 | Threads for ranger and the per-marker step |
| `min.bucket` | 15 | Minimum terminal node size |
| `mtry` | `floor(sqrt(ncol(x)))` | Variables tried per split |
| `imp.blocks` | 1 | Number of independent forests whose importances are summed. Each forest imputes missing values separately. |
| `max.trees` / `min.trees` | `Inf` / 0 | Upper and lower bounds on the trees used per marker |
| `pick.threshold` | 1 | P-value threshold. Only test markers chosen for splits significantly often (1 = test all, e.g. 0.01 = only test markers that were used much more often than expected by chance) |
| `pop.str` | `NULL` | Covariates for population structure |

The trait is standardized internally, and samples with a missing trait are dropped. Missing marker values are filled by sampling from observed values. If fewer than `2 * min.bucket` samples remain, all p-values are returned as 1.
