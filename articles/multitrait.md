# Multi-trait fine-mapping

![](work-in-progress.png)

This is a short tutorial for multi-trait fine-mapping, updated in August
2026. Please note that the multi-trait implementation is under constant
development – moreover, the current implementation is tested for only
two traits, and other number of traits may return an error.

Please report any unexpected behaviour to: ville dot karhunen at
mrc-bsu.cam.ac.uk.

## Package installation and load

Make sure that you are using the most recent version (0.3.0) of the
package.

``` r


# install package
# if needed, install also 'remotes' package via 'install.packages("remotes")'
remotes::install_github("vkarhune/finimom")

# load required package
library(finimom)
```

## Example data

The example data are summary statistics for two simulated traits for 842
variants within 50,000 individuals. Both traits have two causal
variants, one of which is shared with both traits.

``` r


load(url("https://github.com/vkarhune/finimomData/raw/main/multi-finimom_example.RData"))

ls()
#> [1] "causals_y1"      "causals_y2"      "LDmat"           "mafs"           
#> [5] "summarystats_y1" "summarystats_y2"

# summary statistics
head(summarystats_y1)
#>              beta          se
#> [1,] -0.012283147 0.010012730
#> [2,] -0.022047978 0.011903286
#> [3,] -0.007153627 0.010388402
#> [4,]  0.002569805 0.007193564
#> [5,]  0.002972988 0.016465898
#> [6,] -0.012660130 0.009985195
head(summarystats_y2)
#>               beta          se
#> [1,]  0.0061683619 0.010012843
#> [2,] -0.0478840719 0.011901768
#> [3,] -0.0495203583 0.010386090
#> [4,]  0.0005496849 0.007193573
#> [5,] -0.0352351055 0.016465150
#> [6,]  0.0067814593 0.009985310

# causal variants
causals_y1
#> [1] 201 318
causals_y2
#> [1] 201 590
```

The linkage disequilibrium (LD) reference used is an out-of-sample
reference, calculated based on 10,000 individuals. The allele
frequencies for the variants are provided in object ‘mafs’.

``` r


# LD matrix
str(LDmat)
#>  num [1:842, 1:842] 1 -0.0993 -0.1146 -0.2146 -0.0627 ...

# allele frequencies
str(mafs)
#>  num [1:842] 0.1128 0.0762 0.1031 0.2632 0.0385 ...
```

## Fine-mapping

The input data for multi-trait fine-mapping are given as lists:

``` r



# betas
beta <- list(summarystats_y1[,"beta"], summarystats_y2[,"beta"])

# standard errors
se <- list(summarystats_y1[,"se"], summarystats_y2[,"se"])
```

Multi-trait fine-mapping can be run as follows:

``` r


res <- multi_finimom(beta = beta, se = se, eaf = mafs, R = LDmat,
                    n = rep(50000, 2), k = 2, omega = NULL,
                    niter = 62500, burnin = 12500,
                    standardize = TRUE,
                    verbose = TRUE,
                    insampleLD = FALSE,
                    clump = TRUE, clump_r2 = 0.995^2, check_ld = TRUE,
                    u = 2,
                    purity = 0.5,
                    h2cap = NULL, cprior = 0.1, zeta = 0.5)
#> Clumping variants at r2=0.99
#> Sampling from the posterior...
#> 
#> 62500 iterations done in 3.75 seconds

res$sets
#> [[1]]
#> [[1]][[1]]
#> [1]  84 201 288
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 590 598 603 645 653 671
#> 
#> [[2]][[2]]
#> [1]  84 201 288
```

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
