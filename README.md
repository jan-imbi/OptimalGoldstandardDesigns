
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OptimalGoldstandardDesigns <a href='https://github.com/jan-imbi/OptimalGoldstandardDesigns'><img src='man/figures/sticker.png' align="right" height="550" /></a>

<!-- badges: start -->
<!-- badges: end -->

This package implements methods for optimization of three-arm
gold-standard non-inferiority trial design parameters as described in
**TODO: add citation when it appears**

## Installation

The recommended to access the code in this repository is by cloning the
repository via git (download and install
[git](https://git-scm.com/download), open a git shell and type
`git clone https://github.com/jan-imbi/OptimalGoldstandardDesigns`).

You may then open the project in
[Rstudio](https://www.rstudio.com/products/rstudio/download/#download)
and install it by clicking on “Install and Restart” in the “Build” menu.

If you do not wish to review the code and just want to use the functions
within this package, you may install it via the following command:

``` r
remotes::install_github("jan-imbi/OptimalGoldstandardDesigns")
```

## General useage instructions for the functions provided in this package

After installing the package, you made load its function into the
namespace like this:

``` r
# library(OptimalGoldstandardDesigns)
```

**TODO: rewrite examples.**

## Sanity checks and comparison to past works

In the folder *tests*, you can find a subfolder *testthat*. Within this
folder, there are two files. One file is called
*schloemer_diss_appendix.R*, it contains some functions from the
appendix of the dissertation of [Patrick
Schlömer](https://d-nb.info/1072225700/34) (Schlömer 2014).

The other file, *test_sanity_checks.R*, contains some sanity checks for
the presented calculations. For example, it shows that the optimal
versions of designs 2 calculated using our methods agree with the
optimal designs calculated from the functions of Schlömer’s
dissertation. Other checks availble in this file include verification of
the covariance formula and verification that type I errors are properly
controlled. There is also an example which highlights that type I error
is not controlled when using an improper testing procedure with binding
futility boundaries.

## Disclaimer

This package is still under active development and may change at any
time.

At the moment, using this package is not very safe. It is easily
possible to reach an undefined state within the procedures contained in
this package.

This package is not validated for use in clinical trials. If you use
this you use this package to aid in planning a real clinical trial, use
extreme caution. You are responsible to verify the validity of the
results and the sensibility of the trial design. The author of this
package offers absolutely no warranty for any result obtained using this
package. You are using this completly at your own risk.

Feel free to message me if you need assistance in using this package and
in planning a gold-standard non-inferiority trial with optimal design
parameters.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-schlomer2014group" class="csl-entry">

Schlömer, Patrick. 2014. “Group Sequential and Adaptive Designs for
Three-Arm’gold Standard’non-Inferiority Trials.” PhD thesis, Universität
Bremen. <https://d-nb.info/1072225700/34>.

</div>

</div>
