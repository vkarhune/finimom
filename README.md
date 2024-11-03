
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- remember to knit this separately with devtools::build_readme() -->

# finimom <img src="man/figures/finimom.png" align="right" width="125"/>

<!-- badges: start -->
<!-- badges: end -->

> ***NEW: MULTI-TRAIT FINE-MAPPING*** – see
> [tutorial](https://vkarhune.github.io/finimom/articles/multitrait.html).

This R package implements a method for fine-mapping summarised genetic
data using product inverse-moment prior for the effect sizes, and a
beta-binomial prior for the model dimension.

The method is described in:

Karhunen V, Launonen I, Järvelin MR, Sebert S, Sillanpää MJ. Genetic
fine-mapping from summary data using a nonlocal prior improves detection
of multiple causal variants. *Bioinformatics*.
[doi:10.1093/bioinformatics/btad396](https://doi.org/10.1093/bioinformatics/btad396).

## Installation

You can install `finimom` from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("vkarhune/finimom")
```

## Minimal example

Here is a minimal example on using finimom:

``` r

# load package
library(finimom)

# example data:
str(exampledata)
#> List of 9
#>  $ betahat   : num [1:363] 0.0804 0.0927 -0.0253 -0.0165 -0.0128 ...
#>  $ sebetahat : num [1:363] 0.0601 0.0595 0.0968 0.028 0.0276 ...
#>  $ allelefreq: num [1:363] 0.0527 0.0541 0.0194 0.3751 0.4224 ...
#>  $ truebetas : num [1:363] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ causals   : num [1:3] 95 314 329
#>  $ n         : int 2763
#>  $ n_ref     : int 2637
#>  $ insampleLD: num [1:363, 1:363] 1 0.971 0.545 0.288 0.244 ...
#>  $ refLD     : num [1:363, 1:363] 1 0.972 0.558 0.3 0.268 ...

# effect size estimates and their standard errors
beta <- exampledata$betahat
se <- exampledata$sebetahat

# allele frequencies
eaf <- exampledata$allelefreq

# LD matrix
R <- exampledata$insampleLD

# run finimom
res <- finimom(beta = beta, se = se, eaf = eaf, R = R,
               n = exampledata$n,
               insampleLD = TRUE,
               verbose = TRUE, ala = TRUE) # set ala = TRUE for approximate Laplace
#> Calculating tau based on the sample size.
#> Clumping variants at r2=0.98
#> Sampling from the posterior...
#> 
#> 12500 iterations done in 0.97 seconds

# output: credible sets
res$sets
#> [[1]]
#>  [1] 295 300 301 309 310 314 315 316 317 322 327 330 342 345 351 353
#> 
#> [[2]]
#>  [1] 313 320 321 323 325 326 329 331 333 337 339 340 341 343 344 348 350 352 355
#> [20] 360 362
#> 
#> [[3]]
#>  [1]  39  55  63  68  69  72  73  74  75  76  77  78  79  81  82  83  85  86  87
#> [20]  88  90  92  93  94  95  96  97  99 114 116 117 118 124 126 128 129 131 133
#> [39] 136 147 149 150 151 159

# compare with the true causal variants
exampledata$causals
#> [1]  95 314 329

# check whether the true causals were captured in the credible sets
lapply(res$sets, function(x) exampledata$causals %in% x)
#> [[1]]
#> [1] FALSE  TRUE FALSE
#> 
#> [[2]]
#> [1] FALSE FALSE  TRUE
#> 
#> [[3]]
#> [1]  TRUE FALSE FALSE
```

## Further details

Further examples are given in the articles. Please report any bugs or
unexpected behaviour to <ville.karhunen@mrc-bsu.cam.ac.uk>.
