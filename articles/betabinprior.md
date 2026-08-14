# Beta-binomial prior for model dimension

Finimom employs a beta-binomial prior for model dimension $`d`$:

``` math
\mathbb{P}(d = k | a, b) = \binom{p}{k}\frac{B(a + k, p - k + b)}{B(a, b)}, \quad a, b > 0, \quad d = 1, \dots, K,
```

where $`p`$ is the number of variants and $`K`$ the maximum model size.
The priors of this form with $`a = 1`$ and $`b = p^u`$ with $`u > 1`$,
are discussed in Castillo and van der Vaart (2012) and in Castillo et
al. (2015). The parameter $`u`$ here controls the amount of prior
density for smaller models, with larger values of $`u`$ giving more
prior mass to smaller models.

Using a linkage disequilibrium (LD) matrix from an external dataset tend
to increase the false positive rate. In our formulation, parameter $`u`$
provides a flexible way to adjust for this. The default values are
$`u = 2`$ when using in-sample LD matrix, and $`u = 2.25`$ when using
out-of-sample LD matrix.

We demonstrate the prior for model dimension using the example dataset:

``` r

library(finimom)

(p <- length(exampledata$betahat))
#> [1] 363

maxsize <- 10

a <- 1
u <- 1.5

val <- exp(sapply(seq_len(maxsize), dbb, p = p, a = a, b = p^u))
(val <- val/sum(val))
#>  [1] 9.502616e-01 4.727099e-02 2.345333e-03 1.160564e-04 5.727772e-06
#>  [6] 2.819359e-07 1.384076e-08 6.776590e-10 3.309028e-11 1.611478e-12

plot(val, type = "b", ylim = c(0, 1))
```

![](betabinprior_files/figure-html/setup-1.png)

And for different values of $`u`$:

``` r


us <- c(1.05, 1.5, 2, 2.25)

vals <- lapply(us, function(u){
  b <- p^u
  out <- exp(sapply(1:10, dbb, a = a, p = p, b = b))
  out <- out/sum(out)
})

plot(vals[[1]], type = "b", ylim = c(0, 1))
invisible(lapply(2:4, function(i) lines(vals[[i]], type = "b", lty = i)))
```

![](betabinprior_files/figure-html/compareu-1.png)

The same on a log scale:

``` r


plot(vals[[1]], type = "b", log = "y", ylim = range(unlist(vals)))
invisible(lapply(2:4, function(i) lines(vals[[i]], type = "b", lty = i)))
```

![](betabinprior_files/figure-html/compareulog-1.png)

## References

Castillo and van der Vaart (2012). Needles and Straw in a Haystack:
Posterior concentration for possibly sparse sequences. *The Annals of
Statistics*.

Castillo et al. (2015). Bayesian linear regression with sparse priors.
*The Annals of Statistics*.

## Session information

``` r


sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] finimom_0.3.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    Rcpp_1.1.2        yaml_2.3.12       otel_0.2.0       
#> [25] jsonlite_2.0.0    rlang_1.3.0       fs_2.1.0
```
