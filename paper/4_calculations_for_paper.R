library(mvtnorm)
library(nloptr)
library(doFuture)
library(foreach)
library(doRNG)
library(tidyverse)
library(here)
library(future)
source(here("R", "1_conditional_probability_functions.R"))
source(here("R", "2_design_functions.R"))
source(here("R", "3_optimization_methods.R"))

plan(list(
  tweak(multisession, workers = 4) # consider using a sequential plan for low-core machines.
))
set.seed(19102021)

sigma_all <- 1
alpha_all <- .025
maxpts_all <- 4097    # consider using smaller maxpts (used for calculating mvt normal intergrals) for testing
tol_all <- 5e-10      # consider using lower tolerance for testing
maxeval_all <- 20000   # consider using smaller maxeval (used in the optimization routine) for testing
# tol_all <- 1e-1
# maxeval_all <- 2


boundary_conditions_table_1 <- expand.grid(
  type_I_error = c(0.025),
  type_II_error = c(.2, .1),
  Delta = c(0.3),
  alternative_TP = c(0.6),
  lambda = c(1),
  kappa = c(0),
  sigma_T = c(1),
  sigma_P = c(1),
  sigma_C = c(1),
  tol = c(tol_all),
  maxpts = c(maxpts_all),
  maxeval = c(maxeval_all)
) %>% apply(1, as.list)
# ncores = NULL means use as much cores as there are scenarios to simulate.
# Consider limiting ncores on low-core machines.
tab1 <-  future(init_successive_optimum_designs(boundary_conditions_table_1, NULL, 1:6), seed = TRUE)

boundary_conditions_table_2 <- expand.grid(
  type_I_error = c(0.025),
  type_II_error = c(.2, .1),
  Delta = c(0.3),
  alternative_TP = c(0.6),
  lambda = c(1, .9, .8, .7, .6, .5, .4, .3, .2, .1),
  kappa = c(0),
  sigma_T = c(1),
  sigma_P = c(1),
  sigma_C = c(1),
  tol = c(tol_all),
  maxpts = c(maxpts_all),
  maxeval = c(maxeval_all)
) %>% apply(1, as.list)
# ncores = NULL means use as much cores as there are scenarios to simulate.
# Consider limiting ncores on low-core machines.
tab2 <-  future(init_successive_optimum_designs(boundary_conditions_table_2, NULL, c(1, 2, 6)), seed = TRUE)

boundary_conditions_table_3 <- expand.grid(
  type_I_error = c(0.025),
  type_II_error = c(.2, .1),
  Delta = c(0.3),
  alternative_TP = c(0.6),
  lambda = c(1),
  kappa = c(0, .5, 1, 1.5, 2, 2.5, 3),
  sigma_T = c(1),
  sigma_P = c(1),
  sigma_C = c(1),
  tol = c(tol_all),
  maxpts = c(maxpts_all),
  maxeval = c(maxeval_all)
) %>% apply(1, as.list)
# ncores = NULL means use as much cores as there are scenarios to simulate.
# Consider limiting ncores on low-core machines.
tab3 <-  future(init_successive_optimum_designs(boundary_conditions_table_3, NULL, c(1, 2, 6)), seed = TRUE)

boundary_conditions_table_4 <- expand.grid(
  type_I_error = c(0.025),
  type_II_error = c(.2),
  Delta = c(0.3),
  alternative_TP = c(0.6),
  lambda = c(1, .9, .8, .7, .6, .5),
  kappa = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
  sigma_T = c(1),
  sigma_P = c(1),
  sigma_C = c(1),
  tol = c(tol_all),
  maxpts = c(maxpts_all),
  maxeval = c(maxeval_all)
) %>% apply(1, as.list)
# ncores = NULL means use as much cores as there are scenarios to simulate.
# Consider limiting ncores on low-core machines.
tab4 <- future(init_successive_optimum_designs(boundary_conditions_table_4, NULL, c(1, 2, 6)), seed = TRUE)


saveRDS(value(tab1), here("dat", "first_tab.rds"))
saveRDS(value(tab2), here("dat", "second_tab.rds"))
saveRDS(value(tab3), here("dat", "third_tab.rds"))
saveRDS(value(tab4), here("dat", "fourth_tab.rds"))
