# Samples from the posterior distribution

Samples from the posterior distribution

## Usage

``` r
posterior_samples_multitrait(
  beta,
  se,
  eaf,
  R,
  k,
  omega,
  u,
  n,
  maxsize,
  tau0,
  r0,
  niter,
  burnin,
  p,
  seed = 456,
  excl.burnin = TRUE,
  a0 = 1,
  b0 = NULL,
  inds0 = NULL,
  standardize = TRUE,
  verbose = TRUE,
  clump = TRUE,
  clump_r2 = 0.99^2,
  check_ld = FALSE,
  ala = NULL,
  h2cap = NULL,
  lam = cprior,
  collinear = collinear,
  zeta = zeta
)
```

## Arguments

- beta:

  Vector of effect sizes.

- se:

  Vector of standard errors.

- eaf:

  Vector of effect allele frequencies.

- R:

  LD matrix.

- maxsize:

  The maximum number of causal variants.

- tau0:

  Parameter tau.

- r0:

  Parameter \$r\$.

- niter:

  Number of iterations.

- burnin:

  Number of burn-in samples.

- p:

  Number of variants.

- seed:

  Random seed.

- excl.burnin:

  Should the burn-in be excluded?

- a0:

  Hyperparameter a for the model size prior.

- b0:

  Hyperparameter b for the model size prior.

- inds0:

  Initial model indices (not used).

- standardize:

  Should the effect sizes be standardised? Defaults to TRUE.

- verbose:

  Verbose output.

- clump:

  Whether to clump extremely highly correlated variants.

- clump_r2:

  Clumping threshold for extremely highly correlated variants.

- check_ld:

  Should the Z-scores be checked for LD discrepancy?

- ala:

  Whether Approximate Laplace should be used?

## Value

List.
