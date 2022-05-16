library(here)
library(future)
library(testthat)
library(fpCompare)

library(mvtnorm)
library(nloptr)
library(doFuture)
library(foreach)
library(doRNG)

source(here("R", "1_conditional_probability_functions.R"))
source(here("R", "2_design_functions.R"))
source(here("R", "3_optimization_methods.R"))
source(here("R", "table_format_functions.R"))
source(here("tests", "testthat", "schloemer_diss_appendix.R"))

# Read all data produced for paper
first_tab <- readRDS(here("dat", "first_tab.rds"))
second_tab <- readRDS(here("dat", "second_tab.rds"))
third_tab <- readRDS(here("dat", "third_tab.rds"))
fourth_tab <- readRDS(here("dat", "fourth_tab.rds"))

# Create lists containing the optimal designs

## Names of the design columns in the tables
dlist_names <- c(
  "optD_single_stage",
  "optD_no_futility_fixed_c",
  "optD_no_futility",
  "optD_nonbinding_futility",
  "optD_fully_sequential",
  "optD_closed_testing"
)

tab_list <- list(first_tab, second_tab, third_tab, fourth_tab)

design6_list <- lapply(
  tab_list,
  function(x) as.list(x$optD_closed_testing)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_closed_testing"
    return(x)
  })
design5_list <- lapply(
  list(first_tab),
  function(x) as.list(x$optD_fully_sequential)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_fully_sequential"
    return(x)
  })

design4_list <- lapply(
  list(first_tab),
  function(x) as.list(x$optD_nonbinding_futility)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_nonbinding_futility"
    return(x)
  })

design3_list <- lapply(
  list(first_tab),
  function(x) as.list(x$optD_no_futility)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_no_futility"
    return(x)
  })

design2_list <- lapply(
  tab_list,
  function(x) as.list(x$optD_no_futility_fixed_c)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_no_futility_fixed_c"
    return(x)
  })

design1_list <- lapply(
  tab_list,
  function(x) as.list(x$optD_single_stage)
) %>%
  unlist(recursive = FALSE) %>%
  lapply(function(x) {
    x$label <- "optD_single_stage"
    return(x)
  })

all_design_list <- c(
  design1_list,
  design2_list,
  design3_list,
  design4_list,
  design5_list,
  design6_list
)

two_stage_design_list <- c(
  design2_list,
  design3_list,
  design4_list,
  design5_list,
  design6_list
)



# Test that state probabilities of all designs sum to 1
sum_probability_of_states <- function(D, hypothesis = "H0") {
  if (isTRUE(D$nonsequential_futility)) {
    sum(unlist(D$finalStateProbs[[hypothesis]])) -
      # The two probabilities below were calculated for presentational purposes and
      # are not actually state probabilites
      D$finalStateProbs[[hypothesis]]$TP1F -
      D$finalStateProbs[[hypothesis]]$TP1_TC1F
  } else {
    sum(unlist(D$finalStateProbs[[hypothesis]]))
  }
}

test_that("state probabilities of all designs sum to 1", {
  H0_probs <- sapply(all_design_list, sum_probability_of_states, hypothesis = "H0")
  expect_equal(unname(H0_probs), rep(1, length(H0_probs)), tolerance = 1e-6)

  H1_probs <- sapply(all_design_list, sum_probability_of_states, hypothesis = "H1")
  expect_equal(unname(H1_probs), rep(1, length(H0_probs)), tolerance = 1e-6)
})

make_table_from_dlist <- function(dlist) {
  tlist <- list()
  for (di in seq_along(dlist)) {
    d <- dlist[[di]]
    singlestage <- (length(unlist(d$n)) == 3)
    ccc_ <- ccc_wrt_nmax(if (singlestage) {
      d$cc
    } else {
      d$ccc
    }, d$maxn, d$n, singlestage)
    tlist[[length(tlist) + 1]] <- tibble(
      beta = d$type_II_error,
      Delta = d$Delta,
      `Design` = switch(d$label,
        "optD_single_stage" = "1",
        "optD_no_futility_fixed_c" = "2",
        "optD_no_futility" = "3",
        "optD_nonbinding_futility" = "4",
        "optD_fully_sequential" = "5",
        "optD_closed_testing" = "6"
      ),
      # `$c_{1, T}$` =  ccc_[[1]][["T"]],
      # `$c_{1, P}$` =  ccc_[[1]][["P"]],
      # `$c_{1, C}$` =  ccc_[[1]][["C"]],
      # `$c_{2, T}$` =  ccc_[[2]][["T"]],
      # `$c_{2, P}$` =  ccc_[[2]][["P"]],
      # `$c_{2, C}$` =  ccc_[[2]][["C"]],
      # `$n_{max}$` = d$maxn,
      `$n_{1, T}$` = if (singlestage) {
        d$n[[1]][["T"]]
      } else {
        d$cumn[[1]][["T"]]
      },
      `$n_{1, P}$` = if (singlestage) {
        d$n[[1]][["P"]]
      } else {
        d$cumn[[1]][["P"]]
      },
      `$n_{1, C}$` = if (singlestage) {
        d$n[[1]][["C"]]
      } else {
        d$cumn[[1]][["C"]]
      },
      `$n_{2, T}$` = if (singlestage) {
        NA_real_
      } else {
        d$cumn[[2]][["T"]]
      },
      `$n_{2, P}$` = if (singlestage) {
        NA_real_
      } else {
        d$cumn[[2]][["P"]]
      },
      `$n_{2, C}$` = if (singlestage) {
        NA_real_
      } else {
        d$cumn[[2]][["C"]]
      },
      `$N^{}_{H_1}$` = d$ASN$H1,
      `$CP_{\\min}$` = d$min_conditional_power,
      `$b_{1, TP, f}$` = case_when(
        TRUE ~ ifelse(singlestage, NA_real_, round(d$b[[1]][["TP"]][["futility"]], digits = 2))
      ),
      `$b_{1, TP, e}$` = d$b[[1]][["TP"]][["efficacy"]],
      `$b_{2, TP, e}$` = ifelse(singlestage, NA_real_, d$b[[2]][["TP"]][["efficacy"]]),
      `$b_{1, TC, f}$` = case_when(
        TRUE ~ ifelse(singlestage, NA_real_, round(d$b[[1]][["TC"]][["futility"]], digits = 2))
      ),
      `$b_{1, TC, e}$` = d$b[[1]][["TC"]][["efficacy"]],
      `$b_{2, TC, e}$` = ifelse(singlestage, NA_real_, d$b[[2]][["TC"]][["efficacy"]])
    )
  }
  return(bind_rows(tlist))
}

# Design 6 for kappa = 0 and lambda = 1 was calculated multiple times for different tables
# Check whether the results are reasonably similar
D6_lambda1_kappa0_beta02 <- lapply(
  design6_list,
  function(D) {
    if ((D$lambda %==% 1) &&
      (D$kappa %==% 0) &&
      (D$type_II_error %==% 0.2)) {
      return(D)
    } else {
      return(NULL)
    }
  }
)
D6_lambda1_kappa0_beta02[sapply(D6_lambda1_kappa0_beta02, is.null)] <- NULL

tab_duplicate_design6 <- make_table_from_dlist(D6_lambda1_kappa0_beta02)
tab_duplicate_design6_differences <- tab_duplicate_design6 %>%
  summarise(across(-c(Design),
    .fns = ~ max(.x) - min(.x), .names = ".{col}_maxdiff"
  ))

# It's not perfect, but reasonably close
print(tab_duplicate_design6, width = 900)
print(tab_duplicate_design6_differences, width = 900)


# Same goes for design 2
D2_lambda1_kappa0_beta02 <- lapply(
  design2_list,
  function(D) {
    if ((D$lambda %==% 1) &&
      (D$kappa %==% 0) &&
      (D$type_II_error %==% 0.2)) {
      return(D)
    } else {
      return(NULL)
    }
  }
)
D2_lambda1_kappa0_beta02[sapply(D2_lambda1_kappa0_beta02, is.null)] <- NULL



tab_duplicate_design2 <- make_table_from_dlist(D2_lambda1_kappa0_beta02)
tab_duplicate_design2_differences <- tab_duplicate_design2 %>%
  summarise(across(-c(Design),
    .fns = ~ max(.x) - min(.x), .names = ".{col}_maxdiff"
  ))

# Because there are much less parameters to optimize, optimal solutions are
# notably closer to each other in design 2
print(tab_duplicate_design2, width = 900)
print(tab_duplicate_design2_differences, width = 900)



# Check if calculation of Covariance matrix is correct using simulation
# Note that, for inexplicable reasons, i consistently
# denoted the variance by "sigma" instead of "sigma^2"...
set.seed(123)

D <- create_Design(
  type_I_error = 0.025,
  type_II_error = .2,
  Delta = .3,
  alternative_TP = .6,
  lambda = 1,
  kappa = 0,
  sigma_T = 5,
  sigma_C = .4,
  sigma_P = .2,
  tol = 1e-7,
  maxpts = 4097,
  maxeval = 1000
)
D$nonsequential_futility <- TRUE

D$n <- list()
D$n[[2]] <- list()
D$n[[1]][["T"]] <- 123
D$n[[1]][["P"]] <- 23
D$n[[1]][["C"]] <- 321
D$n[[2]][["T"]] <- 432
D$n[[2]][["P"]] <- 42
D$n[[2]][["C"]] <- 432

# Calculate covariance matrix Sigma
D$cc <- calc_c(D$n[[1]][["T"]], D)
D$ccc <- calc_cumc(D)
D$cumn <- calc_cumn(D)
D$rho <- calc_rho(D)
D$Sigma <- calc_Sigma(D)

nsim <- 10000
sample_Z_TP1 <- numeric(nsim)
sample_Z_TP2 <- numeric(nsim)
sample_Z_TC1 <- numeric(nsim)
sample_Z_TC2 <- numeric(nsim)

for (i in seq_len(nsim)) {
  sample_T1 <- rnorm(D$n[[1]][["T"]], mean = D$alternative_TP, sd = sqrt(D$sigma$T))
  sample_T2 <- c(sample_T1, rnorm(D$n[[2]][["T"]], mean = D$alternative_TP, sd = sqrt(D$sigma$T)))
  sample_P1 <- rnorm(D$n[[1]][["P"]], mean = 0, sd = sqrt(D$sigma$P))
  sample_P2 <- c(sample_P1, rnorm(D$n[[2]][["P"]], mean = 0, sd = sqrt(D$sigma$P)))
  sample_C1 <- rnorm(D$n[[1]][["C"]], mean = D$alternative_TP, sd = sqrt(D$sigma$C))
  sample_C2 <- c(sample_C1, rnorm(D$n[[2]][["C"]], mean = D$alternative_TP, sd = sqrt(D$sigma$C)))

  sample_rho_TP1 <- 1 / sqrt(D$sigma$T / D$n[[1]]$T + D$sigma$P / D$n[[1]]$P)

  sample_Z_TP1[i] <- (mean(sample_T1) - mean(sample_P1)) * 1 / sqrt(D$sigma$T / D$n[[1]]$T + D$sigma$P / D$n[[1]]$P)
  sample_Z_TP2[i] <- (mean(sample_T2) - mean(sample_P2)) * 1 / sqrt(D$sigma$T / D$cumn[[2]]$T + D$sigma$P / D$cumn[[2]]$P)
  sample_Z_TC1[i] <- (mean(sample_T1) - mean(sample_C1) + D$Delta) * 1 / sqrt(D$sigma$T / D$n[[1]]$T + D$sigma$C / D$n[[1]]$C)
  sample_Z_TC2[i] <- (mean(sample_T2) - mean(sample_C2) + D$Delta) * 1 / sqrt(D$sigma$T / D$cumn[[2]]$T + D$sigma$C / D$cumn[[2]]$C)
}

cov(data.frame(sample_Z_TP1, sample_Z_TP2, sample_Z_TC1, sample_Z_TC2))
D$Sigma

test_that("Difference between sample variance and calculated variance for n=10000 repetitions is
smaller than <.01", {
  expect_true(all((cov(data.frame(sample_Z_TP1, sample_Z_TP2, sample_Z_TC1, sample_Z_TC2)) -
    D$Sigma) < .01))
})


# Check if type I error is controlled

calc_local_type_I_error_TP <- function(D) {
  # A_ is a list of projection matrices
  A_ <- D$A_
  mu_ <- D$mu_vec[["H0"]]
  Sigma <- D$Sigma
  b <- D$b

  # Designs with no or nonbinding futility boundaries
  if ((b[[1]][["TP"]][["futility"]] == -Inf) || (D$label %in% dlist_names[1:4])) {
    b[[1]][["TP"]][["futility"]] <- -10
  }

  P <- list()
  P[["TP1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1"]] %*% mu_),
    sigma = A_[["TP1"]] %*% Sigma %*% t(A_[["TP1"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]]),
    upper = c(10),
    algorithm = Miwa(steps = D$maxpts)
  )[1]

  P[["TP12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP12"]] %*% mu_),
    sigma = A_[["TP12"]] %*% Sigma %*% t(A_[["TP12"]]),
    lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]]),
    upper = c(b[[1]][["TP"]][["efficacy"]], 10),
    algorithm = Miwa(steps = D$maxpts)
  )[1]
  sum(unlist(P))
}

calc_local_type_I_error_TC <- function(D) {
  # A_ is a list of projection matrices
  A_ <- D$A_
  mu_ <- D$mu_vec[["H0"]]
  Sigma <- D$Sigma
  b <- D$b

  # Designs with no or nonbinding futility boundaries
  if ((b[[1]][["TC"]][["futility"]] == -Inf) || (D$label %in% dlist_names[1:4])) {
    b[[1]][["TC"]][["futility"]] <- -10
  }

  P <- list()
  P[["TC1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TC1"]] %*% mu_),
    sigma = A_[["TC1"]] %*% Sigma %*% t(A_[["TC1"]]),
    lower = c(b[[1]][["TC"]][["efficacy"]]),
    upper = c(10),
    algorithm = Miwa(steps = D$maxpts)
  )[1]

  P[["TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TC12"]] %*% mu_),
    sigma = A_[["TC12"]] %*% Sigma %*% t(A_[["TC12"]]),
    lower = c(b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(b[[1]][["TC"]][["efficacy"]], 10),
    algorithm = Miwa(steps = D$maxpts)
  )[1]
  sum(unlist(P))
}

local_type_I_errors_TP <- sapply(two_stage_design_list, calc_local_type_I_error_TP)
local_type_I_errors_TC <- sapply(two_stage_design_list, calc_local_type_I_error_TC)
test_that("Local type I errors are controlled", {
  expect_lte(max(local_type_I_errors_TP), .025 + 1e-10)
  expect_lte(max(local_type_I_errors_TC), .025 + 1e-10)
})

## Now check global type I error for TC testing problem

calc_global_type_I_error_TC_for_design <- function(mu_TP, D, nonsequential_futility) {
  A_ <- D$A_
  mu_ <- D$mu_vec[["H0"]]
  mu_[1:2] <-
    c(
      mu_TP / D$rho[[1]][["TC"]],
      mu_TP / D$rho[[2]][["TC"]]
    )
  Sigma <- D$Sigma
  b <- D$b

  # Hack to get MiWa algorithm to work with infinite boundaries and reasonable accuracy
  pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
  pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
  nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
  nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

  for (i in seq_len(length(b))) {
    for (j in names(b[[i]])) {
      for (k in names(b[[i]][[j]])) {
        if (b[[i]][[j]][[k]] == Inf) {
          b[[i]][[j]][[k]] <- pInf[[i]][[j]]
        }
        if (b[[i]][[j]][[k]] == -Inf) {
          b[[i]][[j]][[k]] <- nInf[[i]][[j]]
        }
      }
    }
  }
  P <- list()
  P[["TP1E_TC1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
    sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
  )[1]

  P[["TP1E_TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC12"]] %*% mu_),
    sigma = A_[["TP1_TC12"]] %*% Sigma %*% t(A_[["TP1_TC12"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["efficacy"]], pInf[[2]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
  )[1]

  if (nonsequential_futility){
  P[["TP12E_TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP12_TC12"]] %*% mu_),
    sigma = A_[["TP12_TC12"]] %*% Sigma %*% t(A_[["TP12_TC12"]]),
    lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TC"]], pInf[[2]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
  )[1]
  } else{
    P[["TP12E_TC12E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC2"]] %*% mu_),
      sigma = A_[["TP12_TC2"]] %*% Sigma %*% t(A_[["TP12_TC2"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
    )[1]
  }

  sum(unlist(P))
}

calc_max_type_I_error_TC <- function(D, nonsequential_futility=TRUE){
worst_type_I <-  optimize(calc_global_type_I_error_TC_for_design,
         interval = c(0, 10), maximum=TRUE,
         D=D, nonsequential_futility = nonsequential_futility)
worst_type_I$objective
}

# See that global type I error is indeed inflated in design 6 if you test according to the
# fully sequential procedure instead of the hierarchical procedure according to the adapted
# closed testing principle

# Since these calculations are pretty slow, we will only look at some examples
design6_sample <- c(1, 10, 45, 80, 90)
design6_inflations <- sapply(design6_list[design6_sample],
                             calc_max_type_I_error_TC, nonsequential_futility=FALSE)
print("Example of type I error inflation in TC testing problem when not testing according to the closed testing principle in design 6:")
max(design6_inflations)

design6_not_inflated <- sapply(design6_list[design6_sample],
                             calc_max_type_I_error_TC, nonsequential_futility=TRUE)
print("Max type I error with closed testing:")
max(design6_not_inflated)

# Check type I error in TC testing problem for design 5
test_that("Type I error in design 5 is not infalted", {
  design5_TC_errors <- sapply(design5_list,
    calc_max_type_I_error_TC,
    nonsequential_futility = FALSE
  )
  for (i in seq_along(design5_TC_errors)) {
    expect_lt(design5_TC_errors[i], design5_list[[i]]$type_I_error)
  }
})



# Since these calculations are pretty slow, we will only look at some examples
D2_lambda1_kappa0_beta02 <- lapply(
  design2_list,
  function(D) {
    if ((D$lambda %==% 1) &&
        (D$kappa %==% 0)) {
      return(D)
    } else {
      return(NULL)
    }
  }
)
D2_lambda1_kappa0_beta02[sapply(D2_lambda1_kappa0_beta02, is.null)] <- NULL


test_that("Optimal Design 2 are reasonably close", {
  for (i in seq_along(D2_lambda1_kappa0_beta02)){
    d <- D2_lambda1_kappa0_beta02[[i]]
    expect_lt(max(abs(compare_design2_with_schloemer(d))), .01)
  }
})

