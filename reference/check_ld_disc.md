# Check discrepancy of the linkage equilibrium and test statistics

Check discrepancy of the linkage equilibrium and test statistics

## Usage

``` r
check_ld_disc(indices, z, Chi2_quantile = NULL, LDm, clump_r2)
```

## Arguments

- indices:

  Variant indices.

- z:

  Test statistic.

- Chi2_quantile:

  Quantile for declaring a mismatch. The default is the median of a
  Chi-squared distribution with one degree of freedom.

- LDm:

  LD matrix.

- clump_r2:

  Clumping threshold for extremely highly correlated variants.

## Value

List of indices within each LD clump, corrected for the LD discrepancy.
