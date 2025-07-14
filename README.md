
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
data("ddk_2011")
trk <- ddk_2011$tracking
high <- 1-ddk_2011$lowstream
ts <- ddk_2011$totalscore
yy <- (ts - mean(ts[trk==0],na.rm=T))/sd(ts[trk==0],na.rm=T)
xx <- ddk_2011$percentile
yc <- yy[trk==1]
xc <- xx[trk==1]
dc <- high[trk==1]

A <- rd.qte(y=yc,x=xc,d=dc,x0=50,z0=NULL,tau=(1:9/10),bdw=20,cov=0,bias=1)
summary(A,alpha=0.1)
#> 
#> 
#>                                  QTE                                   
#> ---------------------------------------------------------------------- 
#>              Bias cor.    Pointwise         Uniform      
#>     Tau         Est.     Robust S.E.    90% Conf. Band  
#>      0.1      -0.104       0.138      -0.437       0.228
#>      0.2      -0.001       0.145      -0.351       0.349
#>      0.3      -0.068       0.157      -0.447       0.312
#>      0.4      -0.074       0.164      -0.470       0.322
#>      0.5      -0.157       0.185      -0.604       0.290
#>      0.6      -0.069       0.220      -0.599       0.461
#>      0.7      -0.020       0.268      -0.667       0.628
#>      0.8      -0.023       0.308      -0.767       0.721
#>      0.9      -0.003       0.267      -0.646       0.641
```
