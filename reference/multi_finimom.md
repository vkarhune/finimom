# Multi-trait fine-mapping using an inverse-moment prior.

Multi-trait fine-mapping using an inverse-moment prior.

## Usage

``` r
multi_finimom(
  beta,
  se,
  eaf,
  R,
  n = NULL,
  k = NULL,
  omega = NULL,
  cs = TRUE,
  cs_level = 0.95,
  pip = FALSE,
  maxsize = 10,
  tau = NULL,
  r = 1,
  niter = 62500,
  burnin = 12500,
  seed = 456,
  excl.burnin = TRUE,
  a0 = 1,
  b0 = NULL,
  u = NULL,
  inds0 = NULL,
  standardize = TRUE,
  verbose = TRUE,
  insampleLD = NULL,
  clump = TRUE,
  clump_r2 = 0.995^2,
  check_ld = FALSE,
  purity = 0.5,
  h2cap = NULL,
  cprior = NULL,
  collinear = sqrt(0.9),
  zeta = NULL,
  Cmatmethod = "tcor"
)
```

## Arguments

- beta:

  Effect size estimates.

- se:

  Standard errors.

- eaf:

  Effect allele frequencies.

- R:

  LD-matrix.

- n:

  GWAS sample size. If provided, then tau is calculated based on n.

- cs:

  Are credible sets returned?

- cs_level:

  Credible set level.

- pip:

  Are posterior inclusion probabilities returned?

- maxsize:

  Maximum model size.

- tau:

  Prior parameter tau.

- r:

  Prior parameter r.

- niter:

  Number of iterations.

- burnin:

  Number of burn-in iterations.

- seed:

  Random seed.

- excl.burnin:

  Should burn-in be excluded?

- a0:

  Hyperparameter a for the model size prior.

- b0:

  Hyperparameter b for the model size prior (but preferably use
  parameter u).

- u:

  Hyperparameter for model size prior. Defaults to 2 for in-sample LD
  matrix and 2.25 for out-of-sample LD matrix.

- inds0:

  Indices for the starting model.

- standardize:

  Should the effect sizes be standardised? Defaults to TRUE.

- verbose:

  Verbose output.

- insampleLD:

  Is in-sample LD used?

- clump:

  Should clumping be done for extremely highly-correlated variants?

- clump_r2:

  Clumping threshold for extremely highly-correlated variants.

- check_ld:

  Should LD discrepancy check be performed?

- purity:

  Credible set purity, defined as the minimum absolute correlation
  between the variants in a credible set.

- ala:

  Whether Approximate Laplace should be used?

## Value

List.
