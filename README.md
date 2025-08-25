
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QTE.RD

<!-- badges: start -->
<!-- badges: end -->

The goal of **QTE.RD** is to provide comprehensive tools for testing,
estimating, and conducting uniform inference on quantile treatment
effects (QTEs) in sharp regression discontinuity (RD) designs. When
treatment effects vary across covariate-groups, QTE.RD facilitates the
estimation, testing, and visualization of heterogeneous effects by
incorporating covariates and applying the robust bias correction methods
developed by Qu, Yoon, and Perron (2024, <doi:10.1162/rest_a_01168>).

The package is available on CRAN and can be loaded by

``` r
library(QTE.RD)
```

## Example

The following example demonstrates how to use the `rd.qte` function from
the **QTE.RD** package, using data from Duflo, Dupas, and Kremer (2011,
AER). It estimates the quantile treatment effects of tracking on student
achievement.

``` r
data(ddk_2011)
yc <- ddk_2011$ts_std[ddk_2011$tracking==1]
xc <- ddk_2011$percentile[ddk_2011$tracking==1]
dc <- ddk_2011$highstream[ddk_2011$tracking==1]

A <- rd.qte(y=yc,x=xc,d=dc,x0=50,z0=NULL,tau=(1:9/10),bdw=20,bias=1)
summary(A,alpha=0.1)
#> 
#> 
#>                                  QTE                                   
#> ---------------------------------------------------------------------- 
#>              Bias cor.    Pointwise         Uniform      
#>     Tau         Est.     Robust S.E.    90% Conf. Band  
#>      0.1      -0.104       0.139      -0.438       0.229
#>      0.2      -0.001       0.141      -0.340       0.337
#>      0.3      -0.068       0.149      -0.425       0.290
#>      0.4      -0.074       0.155      -0.446       0.298
#>      0.5      -0.157       0.177      -0.581       0.268
#>      0.6      -0.069       0.213      -0.581       0.442
#>      0.7      -0.020       0.270      -0.668       0.629
#>      0.8      -0.023       0.318      -0.787       0.740
#>      0.9      -0.003       0.274      -0.661       0.656
```
