# Title

Title

## Usage

``` r
select_causals(
  numcausals,
  LDmat,
  mincorr = 0.5,
  maxcorr = 0.8,
  minanycorr = 0.5
)
```

## Arguments

- numcausals:

  Number of causal variants.

- LDmat:

  Linkage disequilibrium (LD) matrix, ie. correlation matrix of the
  variants.

- mincorr:

  Minimum absolute correlation between a pair of causal variants.

- maxcorr:

  Maximum absolute correlation between a pair of causal variants.

- minanycorr:

  Minimum absolute correlation of a causal variant with at least one
  other variant.

## Value

Indices of causal variants.
