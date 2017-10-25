<!-- README.md is generated from README.Rmd. Please edit that file -->
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
