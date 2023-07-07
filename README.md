
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- remember to knit this separately with devtools::build_readme() -->

# finimom <img src="man/figures/finimom.png" align="right" width="125"/>

<!-- badges: start -->
<!-- badges: end -->

This R package implements a method for fine-mapping summarised genetic
data using product inverse-moment prior for the effect sizes, and a
beta-binomial prior for the model dimension.

The method is described in a pre-print:

Karhunen V, Launonen I, Jarvelin MR, Sebert S, Sillanpaa MJ. Genetic
fine-mapping from summary data using a non-local prior improves
detection of multiple causal variants. *bioRxiv*. 2022:2022-12.
[doi:10.1101/2022.12.02.518898](https://doi.org/10.1101/2022.12.02.518898).

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
#>  $ betahat   : num [1:363] -0.1326 -0.13 -0.1328 0.0096 -0.0422 ...
#>  $ sebetahat : num [1:363] 0.0601 0.0594 0.0967 0.028 0.0276 ...
#>  $ allelefreq: num [1:363] 0.0527 0.0541 0.0194 0.3751 0.4224 ...
#>  $ truebetas : num [1:363] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ causals   : num [1:3] 14 179 195
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
#> 12500 iterations done in 0.72 seconds

# output: credible sets
res$sets
#> [[1]]
#>  [1] 167 169 172 174 177 178 179 180 183 185 186 187 188 190 193 194 197 206 208
#> [20] 209 211 215 223 234 241 245 279
#> 
#> [[2]]
#>  [1] 168 175 178 181 189 195 196 200 201 202 204 205 210 211 212 213 214 215 216
#> [20] 218 219 220 221 223 224 227 229 230 232 234 235 238 239 240 241 242 248 251
#> [39] 254

# compare with the true causal variants
exampledata$causals
#> [1]  14 179 195

# check whether the true causals were captured in the credible sets
lapply(res$sets, function(x) exampledata$causals %in% x)
#> [[1]]
#> [1] FALSE  TRUE FALSE
#> 
#> [[2]]
#> [1] FALSE FALSE  TRUE
```

## Further details

Further examples are given in the articles.
