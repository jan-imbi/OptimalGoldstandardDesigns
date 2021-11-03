
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OptimalGoldstandardDesigns

<!-- badges: start -->
<!-- badges: end -->

This package contains the code used in the calculations for our
manuscript on the optimization of the two-stage group sequential
three-arm gold-standard design for non-inferiority trials.

## Information for Reviewers

In the following, we will provide some pointers which will hopefully
make the reviewing process more convenient.

## Current state of development of this package

You are currently on the *master* branch, the version of the package
presented here reflects the state at submission. Further development of
the package is taking place on the development branch (*dev*) and will
be merged into the master some time in the future.

Goals for future versions are the ability to have more than 2 stages and
an option to use t-test statistics instead of Z-tests, as well as
improvements to documentation and acceleration of the optimization
procedure.

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

## Data used for the preparation of the manuscript

Data used in the generation of the four tables presented in the
manuscript can be found in the *dat* folder.

This data is generated via the script *4_calculations_for_paper.R*,
which is located in the *paper* folder. Execution of this script may
take a long time depending on your machine. If run as it is, the script
will start 80 parallel threads. The server we used for the calculations
was relatively modern as of 2021 and had 50 cores available. On this
machine, the calculation takes about 1 day.

## Table generation and data extraction

The code used for generating the tables and for extracting the numbers
referenced in the manuscript can be found in the *paper* folder. It is
contained in the file *paper_tables_and_numbers.Rmd*.

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

## General useage instructions for the functions provided in this package

### An example featuring ‘Design 6’

After installing the package, you made load its function into the
namespace like this:

``` r
library(OptimalGoldstandardDesigns)
#> Loading required package: nloptr
#> Loading required package: mvtnorm
#> Loading required package: doFuture
#> Loading required package: foreach
#> Loading required package: future
#> Loading required package: doRNG
#> Loading required package: rngtools
```

First, you need to specify some boundary conditions for your
gold-standard design:

``` r
design_template <- create_Design(
  type_I_error = 0.025,
  type_II_error = 0.2,
  Delta = 0.3,
  alternative_TP = 0.6,
  lambda = 1,
  kappa = 0,
  sigma_T = 1,
  sigma_P = 1,
  sigma_C = 1,
  tol = 1e-3,
  maxpts = 500,
  maxeval = 400
  )
```

You can than find optimal design parameters for this design like this:

``` r
print("This is gonna take 5 minutes")
#> [1] "This is gonna take 5 minutes"
optimal_design6_parm <- opt_objective_closed_testing(D = design_template)
```

The output of the function call above is the return value of an `nloptr`
routine. It contains optimized design parameters in the parameterization
we are working in, but not much else. To get some more details, you can
evaluate the objective function at these optimal design parameters:

``` r
optimal_design6 <- objective_closed_testing(x = optimal_design6_parm$solution,
                         D = design_template,
                         return_everything = TRUE)
```

Now we can have a look at some of the characteristics of this optimized
design. Here is the average sample size under the alternative
hypothesis:

``` r
optimal_design6$ASN$H1
#> [1] 339.0018
```

Here are the number of patients to be recruited in the respective stages
and study arms:

``` r
optimal_design6$n
#> [[1]]
#> [[1]]$T
#> [1] 104.5327
#> 
#> [[1]]$P
#> [1] 37.03078
#> 
#> [[1]]$C
#> [1] 107.8289
#> 
#> 
#> [[2]]
#> [[2]]$T
#> [1] 88.77034
#> 
#> [[2]]$P
#> [1] 58.28121
#> 
#> [[2]]$C
#> [1] 91.60705
```

Here are the efficacy and futility boundaries:

``` r
optimal_design6$b
#> [[1]]
#> [[1]]$TP
#> [[1]]$TP$futility
#> [1] 0.2747707
#> 
#> [[1]]$TP$efficacy
#> [1] 2.023543
#> 
#> 
#> [[1]]$TC
#> [[1]]$TC$futility
#> [1] 0.8478128
#> 
#> [[1]]$TC$efficacy
#> [1] 2.253281
#> 
#> 
#> 
#> [[2]]
#> [[2]]$TP
#> [[2]]$TP$efficacy
#> [1] 2.502896
#> 
#> 
#> [[2]]$TC
#> [[2]]$TC$efficacy
#> [1] 2.057282
```

To speed up computation in this tutorial, we chose high tolerances and a
low number of maximum iterations. We can see how this optimum compares
to the designs calculated for the manuscript:

``` r
first_tab <- readRDS(here::here("dat", "first_tab.rds"))
design6_manuscript <- first_tab$optD_closed_testing[[1]]

diff <- optimal_design6$ASN$H1 - design6_manuscript$ASN$H1
print(diff)
#> [1] -0.621878
```

This example turns out to be not that bad. Our calculation with way
smaller tolerances and way more iterations has an advantage in terms of
average sample size of only 0.622 subjects.

We have to remember though that we used relatively high tolerances in
the calculations. Let’s look at the power of this design with a bit more
accuracy:

``` r
check_power_design6 <- optimal_design6
check_power_design6$maxpts <- 4097
check_power_design6$tol <-  1e-8
calc_prob_reject_both("H1", check_power_design6)
#> [1] 0.7997021
```

We can also check out the local type I error for the Treatment vs
Placebo testing problem:

``` r
alphaTP1_design6 <- pmvnorm_(
  mean = as.vector(optimal_design6$A_[["TP1"]] %*% optimal_design6$mu_vec[["H0"]]),
  sigma = optimal_design6$A_[["TP1"]] %*% optimal_design6$Sigma %*% t(optimal_design6$A_[["TP1"]]),
  lower = optimal_design6$b[[1]][["TP"]][["efficacy"]],
  upper = 10,
  algorithm = Miwa(steps = 4097)
)[1]
alphaTP2_design6 <- pmvnorm_(
  mean = as.vector(optimal_design6$A_[["TP12"]] %*% optimal_design6$mu_vec[["H0"]]),
  sigma = optimal_design6$A_[["TP12"]] %*% optimal_design6$Sigma %*% t(optimal_design6$A_[["TP12"]]),
  lower = c(optimal_design6$b[[1]][["TP"]][["futility"]], optimal_design6$b[[2]][["TP"]][["efficacy"]]),
  upper = c(optimal_design6$b[[1]][["TP"]][["efficacy"]], 10),
  algorithm = Miwa(steps = 4097)
)[1]
alphaTP1_design6 + alphaTP2_design6
#>      upper 
#> 0.02501124
```

And the local error for the Treatment vs Control testing problem:

``` r
alphaTC1_design6 <- pmvnorm_(
  mean = as.vector(optimal_design6$A_[["TC1"]] %*% optimal_design6$mu_vec[["H0"]]),
  sigma = optimal_design6$A_[["TC1"]] %*% optimal_design6$Sigma %*% t(optimal_design6$A_[["TC1"]]),
  lower = optimal_design6$b[[1]][["TC"]][["efficacy"]],
  upper = 10,
  algorithm = Miwa(steps = 4097)
)[1]
alphaTC2_design6 <- pmvnorm_(
  mean = as.vector(optimal_design6$A_[["TC12"]] %*% optimal_design6$mu_vec[["H0"]]),
  sigma = optimal_design6$A_[["TC12"]] %*% optimal_design6$Sigma %*% t(optimal_design6$A_[["TC12"]]),
  lower = c(optimal_design6$b[[1]][["TC"]][["futility"]], optimal_design6$b[[2]][["TC"]][["efficacy"]]),
  upper = c(optimal_design6$b[[1]][["TC"]][["efficacy"]], 10),
  algorithm = Miwa(steps = 4097)
)[1]
alphaTC1_design6 + alphaTC2_design6
#>      upper 
#> 0.02501092
```

Calculating power and local type I errors with higher accuracy gives
that this low accuracy example design (very) slightly undercuts the
target power and has slighty inflated local type I errors.

The design above is what is dubbed ‘design 6’ in our manuscript. It is a
‘gold-standard’ three-arm fixed non-inferiority margin design featuring
a group-sequential recruitment process and testing procedure. The
testing procedure includes binding futility boundaries and type I error
recycling and is based on a variation of the closed testing procedure.

### An example featuring ‘Design 4’

Design 4 is similar to Design 6, though futility boundaries are not
binding and no type I error recycling takes place.

``` r
print("This is gonna take 5 minutes")
#> [1] "This is gonna take 5 minutes"
optimal_design4_parm <- opt_objective_nonbinding_futility(D = design_template)
optimal_design4 <- objective_nonbinding_futility(x = optimal_design4_parm$solution,
                         D = design_template,
                         return_everything = TRUE)
```

``` r
optimal_design4$ASN$H1
#> [1] 342.9316
```

``` r
optimal_design4$n
#> [[1]]
#> [[1]]$T
#> [1] 108.9703
#> 
#> [[1]]$P
#> [1] 35.68632
#> 
#> [[1]]$C
#> [1] 107.6134
#> 
#> 
#> [[2]]
#> [[2]]$T
#> [1] 82.90757
#> 
#> [[2]]$P
#> [1] 56.94409
#> 
#> [[2]]$C
#> [1] 81.83024
```

``` r
optimal_design4$b
#> [[1]]
#> [[1]]$TP
#> [[1]]$TP$futility
#> [1] -0.5267136
#> 
#> [[1]]$TP$efficacy
#> [1] 2.027916
#> 
#> 
#> [[1]]$TC
#> [[1]]$TC$futility
#> [1] 0.6018754
#> 
#> [[1]]$TC$efficacy
#> [1] 2.276344
#> 
#> 
#> 
#> [[2]]
#> [[2]]$TP
#> [[2]]$TP$efficacy
#> [1] 2.495381
#> 
#> 
#> [[2]]$TC
#> [[2]]$TC$efficacy
#> [1] 2.084721
```

We can once again compare this to the optimum found in in the
calculations for the paper:

``` r
design4_manuscript <- first_tab$optD_nonbinding_futility[[1]]

diff <- optimal_design4$ASN$H1 - design4_manuscript$ASN$H1
print(diff)
#> [1] -0.8048035
```

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
