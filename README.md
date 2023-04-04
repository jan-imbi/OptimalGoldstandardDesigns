
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OptimalGoldstandardDesigns <a href='https://github.com/jan-imbi/OptimalGoldstandardDesigns'><img src='man/figures/sticker.png' align="right" height="550" /></a>

<!-- badges: start -->

[![doi](https://img.shields.io/badge/doi-10.1002%2Fsim.9630-blue)](https://doi.org/10.1002/sim.9630)
[![Codecov test
coverage](https://codecov.io/gh/jan-imbi/OptimalGoldstandardDesigns/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jan-imbi/OptimalGoldstandardDesigns?branch=master)
[![R-CMD-check](https://github.com/jan-imbi/OptimalGoldstandardDesigns/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jan-imbi/OptimalGoldstandardDesigns/actions/workflows/R-CMD-check.yaml)
[![License](https://img.shields.io/github/license/jan-imbi/OptimalGoldstandardDesigns)](https://github.com/jan-imbi/OptimalGoldstandardDesigns/blob/master/LICENSE.md)
<!-- badges: end -->

This package contains the code used in the calculations for our
[paper](https://doi.org/10.1002/sim.9630) on the optimization of the
two-stage group sequential three-arm gold-standard design for
non-inferiority trials.

It allows for the simultaneous optimization of the allocation ratios for
both stages (in two-stage designs), the efficacy boundaries, and the
futility boundaries. The optimization is performed under type I and II
error constraints and the objective function is customizable by the
user. Methods to optimize two- and one-stage designs are available.

# Installation

You can install the CRAN version of this package by typing the following
command into your R console:

``` r
install.packages("OptimalGoldstandardDesigns")
```

You can install the GitHub Version by typing:

``` r
remotes::install_github("jan-imbi/OptimalGoldstandardDesigns")
```

You can also [clone this
repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
directly from github. This will give you access to the /data/
subdirectory, which contains code to reproduce the examples from the
paper.

# Documentation

You can check out an online version of the documentation here:
[jan-imbi.github.io/OptimalGoldstandardDesigns](https://jan-imbi.github.io/OptimalGoldstandardDesigns/).
You may be interested in checking out the [*Usage guidance*
article](https://jan-imbi.github.io/OptimalGoldstandardDesigns/articles/Introduction.html).

# References

Meis, J, Pilz, M, Herrmann, C, Bokelmann, B, Rauch, G, Kieser, M.
Optimization of the two-stage group sequential three-arm gold-standard
design for non-inferiority trials. *Statistics in Medicine.* 2023; 42(
4): 536– 558. [doi:10.1002/sim.9630](https://doi.org/10.1002/sim.9630).
