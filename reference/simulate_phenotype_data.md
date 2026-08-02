# Title

Title

## Usage

``` r
simulate_phenotype_data(
  X,
  N,
  p,
  causals,
  R2,
  seed,
  meanbeta = 0,
  sdbeta = 1,
  minpower = NULL,
  alpha = NULL
)
```

## Arguments

- X:

  Genotype matrix.

- N:

  Number of individuals.

- p:

  Number of genetic variants.

- causals:

  Indices of causal variants.

- R2:

  Variance explained by causal variants.

- seed:

  Random seed.

- meanbeta:

  Mean of the Gaussian where effect sizes are drawn from.

- sdbeta:

  Standard deviation of the Gaussian where the effect sizes are drawn
  from.

- alpha:

## Value

List of the phenotype and the betas, after scaling to var(Y) = 1.
