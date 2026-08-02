# Extract a subchain from exsisting posterior samples, excluding a burn-in.

Extract a subchain from exsisting posterior samples, excluding a
burn-in.

## Usage

``` r
extract_subchain(posterior_sample, burnin = NULL, niter = NULL)
```

## Arguments

- posterior_sample:

  Object from posterior_samples().

- burnin:

  The number of initial iterations to be excluded.

- niter:

  The number of iterations kept.

## Value

Posterior samples with burn-in excluded.
