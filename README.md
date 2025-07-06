# QTE_RD
Quantile Treatment Effects under the Regression Discontinuity Design

The R package QTE.RD provides comprehensive methods for testing, estimating, and conducting uniform inference on quantile treatment effects (QTEs) in sharp regression discontinuity (RD) designs, incorporating covariates and implementing robust bias correction methods of Qu, Yoon, Perron (2024) doi:10.1162/rest_a_01168.

The package is available on CRAN and can be loaded by

library(QTE.RD)

The package QTE.RD includes four main functions:

rd.qte() estimates QTEs and provides uniform confidence bands, with or without covariates, and with or without robust bias correction.
rdq.test() conducts tests for three hypotheses, related to the significance of treatment effects, homogeneous treatment effects, and uniformly positive or negative treatment effects.
rdq.bandwidth() implements two bandwidth selection rules: the cross-validation bandwidth and the MSE optimal bandwidth.
plot.qte() generates figures summarizing the treatment effects along with their confidence bands.
