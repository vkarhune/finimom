# Obtain credible sets from the posterior samples.

Obtain credible sets from the posterior samples.

## Usage

``` r
get_credible_sets(samples, num_signals, level = 0.95, purity = purity, R = R)
```

## Arguments

- samples:

  Posterior samples from posterior_samples() function.

- num_signals:

  Number of signals.

- level:

  Probability level, default = 0.95.

- purity:

  Credible set purity.

- R:

  Correlation matrix.

## Value

List of credible sets.
