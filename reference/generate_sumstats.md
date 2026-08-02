# Title

Title

## Usage

``` r
generate_sumstats(phenotype, X, adj = NULL)
```

## Arguments

- phenotype:

  Phenotype.

- X:

  Genotype matrix.

- adj:

  Adjustments (defaults to NULL).

## Value

Summary statistics with four columns: Effect size, its standard error,
log(p-value), and p-value (note: min(p) = 1e-300, more accurate given by
log(p-value)).
