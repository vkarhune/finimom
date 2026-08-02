# Estimate tau based on sample size and z score

Estimate tau based on sample size and z score

## Usage

``` r
estimate_tau(n = NULL, q = 0.05, pval = 0.001, z = NULL)
```

## Arguments

- n:

  Sample size in GWAS summary statistics

- q:

  Quantile for the proposed Z score

- pval:

  P-value threshold for detectable effect size (z-score used if both
  p-value and z-score provided)

- z:

  Z-score for detectable effect size (z-score used if both p-value and
  z-score provided)
