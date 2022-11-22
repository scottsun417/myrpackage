
# rpackage625

<!-- badges: start -->
[![R-CMD-check](https://github.com/scottsun417/myrpackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scottsun417/myrpackage/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/scottsun417/myrpackage/branch/main/graph/badge.svg)](https://app.codecov.io/gh/scottsun417/myrpackage?branch=main)
<!-- badges: end -->

The goal of rpackage625 is to implement a logistic regression used for modelling categorical data, by calling an own writing R function named 'manual_logistic_regression'. This function assuming the response variable Y falls into exactly two categories, Y = 1 and the other by Y = 0.

## Installation

To use the package, you can install the development version of rpackage625 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("scottsun417/myrpackage")
```

## Example

This is a basic example which shows you how to use the function manual_logistic_regression in rpackage625:

``` r
# library the rpackage625
library(rpackage625)

# data simulation
set.seed(2022)
x1 = rnorm(30,3,2) + 0.1*c(1:30)
x2 = rbinom(30, 1,0.3)
x3 = rpois(n = 30, lambda = 4)
x3[16:30] = x3[16:30] - rpois(n = 15, lambda = 2)
x = cbind(x1,x2,x3)
y = c(rbinom(5, 1,0.1),rbinom(10, 1,0.25),rbinom(10, 1,0.75),rbinom(5, 1,0.9))

# conduct logistic regression
manual_glm <- manual_logistic_regression(x, y)

# get estimated coefficients
manual_glm$coefficients

# get covariance_matrix of the logistic regression
manual_glm$covariance_matrix

# get Null Deviance
manual_glm$Null_deviance
```

## Getting Help
If you have any questions or would like to discuss further, please feel free to contact me (zhiyisun@umich.edu).

