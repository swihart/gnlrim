<!-- README.md is generated from README.Rmd. Please edit the .Rmd -->
[![Travis-CI Build Status](https://travis-ci.org/swihart/gnlrim.svg?branch=master)](https://travis-ci.org/swihart/gnlrim)

gnlrim
======

The goal of gnlrim is to have a freestanding, smaller footprint instance of `repeated::gnlmix()`. That is, I took the minimum amount of code from `rmutil` and `repeated` to get gnlmix working so that `gnlrim` did not require dependence/imports from `rmutil` and `repeated`. I am going to be making edits on `gnlmix()` to transform it into `gnlrim()` so that it is basically `gnlmix()` with the added functionality:

-   Stand-alone / smaller footprint
-   Use nlminb for constrained optimization
-   Build-in random intercept distributions that are bridging distributions
-   Add a `data=` argument so that we can do away with using `attach()`
-   **Possibly** add the ability to extract random intercepts / do prediction
-   Make some pre-canned, user-friendly calls for marginalized random intercept models (MRIMs)

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

Errors and Their (Potential) Fixes
----------------------------------

       Likelihood returns Inf or NAs: invalid initial values, wrong model, or probabilities too small to calculate

This is sometimes due to not using a named lists for `pmu` and `pmix`.

``` r
## without named lists for `pmu` and `pmix`:
# fit_PPN_err <- gnlrim(y=y_cbind,
#                       mu=~pnorm(a+b*dose+rand),
#                       pmu=c(0,0),
#                       pmix=0.5,
#                       p_uppb = c(Inf ,  Inf,   1-1e-5),
#                       p_lowb = c(-Inf, -Inf,   0+1e-5),
#                       distribution="binomial",
#                       nest=id,
#                       random="rand",
#                       mixture="normal-phi")
# Error in gnlrim(y = y_cbind, mu = ~pnorm(a + b * dose + rand), pmu = c(0,  : 
#   Likelihood returns Inf or NAs: invalid initial values, wrong model, or probabilities too small to calculate

## with named lists for `pmu` and `pmix`:
fit_PPN_fix <- gnlrim(y=y_cbind,
                      mu=~pnorm(a+b*dose+rand),
                      pmu=c(a=0,b=0),
                      pmix=c(phi=0.5),
                      p_uppb = c(Inf ,  Inf,   1-1e-5),
                      p_lowb = c(-Inf, -Inf,   0+1e-5),
                      distribution="binomial",
                      nest=id,
                      random="rand",
                      mixture="normal-phi")
## a, b, phi:
fit_PPN_fix$coeff
#> [1] -0.4832740  0.3560651  0.4203699
```
