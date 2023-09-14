# Running this script will reproduce all results presented within this dissertation.
# The calculations are started on multiple threads and can take advantage of
# multi-core machines. At the time of writing, the calculations took about 4 days to
# complete on a 50-core machine.

# The packages developed as part of this dissertation will likely be extended in the
# future. To ensure reproducibility of the presented results, the current state of the
# packages is archived in the respective github repositories on a branch called
# 'dissertation'. The first step to reproduce the result is to install these specific
# versions of the packages, which is done in the following.
library(remotes)
install_github("jan-imbi/adoptr", ref = "dissertation")
install_github("jan-imbi/adestr", ref = "dissertation")
install_github("jan-imbi/OptimalGoldstandardDesigns", ref = "dissertation")

# Load package the packages developed within this thesis into the environment.
library(adestr)
library(OptimalGoldstandardDesigns)
# Package for data manipulation
library(dplyr)
# The future package provides access to the parallel execution back end.
library(future.apply)
# Set up the parameters for parallel execution. You can lower the
# workers count to limit the number of cores to be used.
plan(multisession, workers = 45L)
# Library to display progress bars and set up of those bars
library(progressr)
handlers(global = TRUE)
handlers("progress")

# The results about to be created will be saved in the "data" sub-directory of the
# current working directory. In case this directory does not exist, it will be
# created now.
dir.create("data")

######################################################################################
##               Point estimates, confidence intervals and p-values                 ##
######################################################################################

# The designs used in the comparison of the presented statistics will be
# derived here. The designs parameters are optimized via the adoptr package.
# The first pair of designs is optimized with assuming a normal endpoint with known
# variance.
datadist <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
# The power is calculated under the point hypothesis mu=0.4.
H_1 <- PointMassPrior(.4, 1)
ess_H0 <- ExpectedSampleSize(datadist, H_0)
ess_H1 <- ExpectedSampleSize(datadist, H_1)
toer  <- Power(datadist, H_0)
power <- Power(datadist, H_1)
# Here, the intial design parameters from which the optimization will start are set.
initial_gs <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "group-sequential",
  dist  = datadist
)
initial_ad <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist
)
# Here, the parameters for the adaptive design are optimized.
optad <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_ad,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 40000)
)
# Caching the spline functions improves performance.
designad <- optad$design
attr(designad, "label") <- "Adaptive"
# Here, the parameters of the group-sequential design are optimized.
optgs <- minimize(
  ess_H1,
  subject_to(
    power >= 0.8,
    toer  <= .025
  ),
  initial_gs,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 40000)
)
designgs <- optgs$design
attr(designgs, "label") <- "Group-sequential"

# Next, the two designs used for the presentation of the statistics in the two-armed
# known variance case are optimized.
datadist_t <- Student(two_armed = TRUE)
ess_H0_t <- ExpectedSampleSize(datadist_t, H_0)
ess_H1_t <- ExpectedSampleSize(datadist_t, H_1)
toer_t  <- Power(datadist_t, H_0)
power_t <- Power(datadist_t, H_1)
initial_gs_t <- get_initial_design(
  theta = .3,
  alpha = .025,
  beta  = .2,
  type_design  = "group-sequential",
  dist  = datadist_t
)
initial_ad_t <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist_t
)
optad_t <- minimize(
  ess_H1_t,
  subject_to(
    power_t >= 0.8,
    toer_t  <= .025
  ),
  initial_ad_t,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 40000)
)
designad_t <- optad_t$design
attr(designad_t, "label") <- "Adaptive"
optgs_t <- minimize(
  ess_H1_t,
  subject_to(
    power_t >= 0.8,
    toer_t  <= .025
  ),
  initial_gs_t,
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-06, maxeval = 40000)
)
designgs_t <- optgs_t$design
attr(designgs_t, "label") <- "Group-sequential"

# Save all of the designs to disk.
saveRDS(designad, "data/designad.rds")
saveRDS(designgs, "data/designgs.rds")
saveRDS(designad_t, "data/designad_t.rds")
saveRDS(designgs_t, "data/designgs_t.rds")

# Here, the parameter lists for the scenarios to be evaluated are defined.
muvec <- sort(c(seq(-0.75, 1.32, .03),
                mean(c(designad@c1f, designad@c1e))/sqrt(designad@n1),
                mean(c(designgs@c1f, designgs@c1e))/sqrt(designgs@n1)))
# The mu grid is saved.
saveRDS(muvec, "data/muvec.rds")
point_estimator_score_list <- list(MSE(), OverestimationProbability())
interval_estimator_score_list <- list(Coverage(), Width(), TestAgreement())
point_estimator_list <- list(SampleMean(),
                             FirstStageSampleMean(),
                             PseudoRaoBlackwell(),
                             BiasReduced(),
                             MinimizePeakVariance(),
                             MedianUnbiasedStagewiseCombinationFunctionOrdering(),
                             MedianUnbiasedMLEOrdering(),
                             MedianUnbiasedLikelihoodRatioOrdering(),
                             MedianUnbiasedScoreTestOrdering())
interval_estimator_list <- list(StagewiseCombinationFunctionOrderingCI(),
                                MLEOrderingCI(),
                                LikelihoodRatioOrderingCI(),
                                ScoreTestOrderingCI(),
                                NeymanPearsonOrderingCI(),
                                RepeatedCI(),
                                NaiveCI())
dist <- list(Normal(two_armed = FALSE))
dist_t <- list(Student(two_armed = TRUE))
designs <- list(designad, designgs)
designs_t <- list(designad_t, designgs_t)
# Set options to control the accuracy of the integration routines.
options(list(
  # Root finding inside estimators
  adestr_tol_roots = 1e-3,
  adestr_maxiter_roots = 1e3,
  # Integrals used inside estimators
  adestr_tol_inner = 5e-3,
  adestr_maxEval_inner = 1e4,
  adestr_absError_inner = 1e-5,
  # Integrals to evaluate estimators
  adestr_tol_outer = 5e-6,
  adestr_maxEval_outer = 3e5,
  adestr_absError_outer = 1e-8
))

### Set seed for just to be save (although basically nothing is random)
seed <- 01062023
set.seed(seed)

# Calculate the integrals of the different scenarios in parallel. This part is where
# most of the computation happens.
result_list <- evaluate_scenarios_parallel(
  score_lists = list(point_estimator_score_list, interval_estimator_score_list,
    point_estimator_score_list, interval_estimator_score_list),
  estimator_lists = list(point_estimator_list, interval_estimator_list,
    point_estimator_list, interval_estimator_list),
  data_distribution_lists = list(dist, dist, dist_t, dist_t),
  design_lists = list(designs, designs, designs_t, designs_t),
  mu_lists = list(muvec, muvec, muvec, muvec),
  sigma_lists = list(1, 1, 1, 1)
)
est_results <- result_list[[1]] |> select(-"EstimatorScoreResult")
ci_results <- result_list[[2]] |> select(-"EstimatorScoreResult")
est_results_t <- result_list[[3]] |> select(-"EstimatorScoreResult")
ci_results_t <- result_list[[4]] |> select(-"EstimatorScoreResult")
# Save the results on disk.
saveRDS(est_results, "data/est_results.rds")
saveRDS(ci_results, "data/ci_results.rds")
saveRDS(est_results_t, "data/est_results_t.rds")
saveRDS(ci_results_t, "data/ci_results_t.rds")

# Clean up computing environment.
rm(datadist, datadist_t, designad, designad_t, designgs, designgs_t, designs,
  designs_t, dist, dist_t, ess_H0, ess_H0_t, ess_H1, ess_H1_t, H_0, H_1, initial_ad,
  initial_ad_t, initial_gs, initial_gs_t, interval_estimator_list,
  interval_estimator_score_list, muvec, optad, optad_t, optgs, optgs_t,
  point_estimator_list, point_estimator_score_list, power, power_t, toer, toer_t)

######################################################################################
##                     Gold-standard non-inferiority designs                        ##
######################################################################################

# Some algorithms, such as MLSL, are non-deterministic. For reproducibility, a seed
# for the random number generator is set.
seed <- 01062023
set.seed(seed)

# Helper function to create outer products of lists of input parameters.
expand.grid.list <- function(...) apply(expand.grid(...), 1, as.list)

# Design parameters for the single-stage design. Note that only parameters which
# deviate from the default parameters need to be supplied, see
# ?optimize_design_onestage for the defaults.
params_onestage <- expand.grid.list(
  nr = 1,
  beta = c(.2, .1),
  print_progress = FALSE,
  round_n = FALSE
)

# Set up progress bar
with_progress({
prog <- progressor(steps = length(params_onestage))
# Optimize the single-stage designs. Note that it is necessary to do this first,
# because the Schloemer Brannath design (design 2) is partially constructed from
# parameters from the optimal single-stage design.
results_onestage <- future_lapply(
  params_onestage,
  function(x) {
    prog("")
    do.call(optimize_design_onestage, x, quote = TRUE)
    },
  future.seed = seed) |>
  OptimalGoldstandardDesigns:::make_table()
})

# These are the design parameters which will be shared across all two-stage designs.
shared_params_twostage <- list(
  round_n = FALSE,
  print_progress = FALSE,
  inner_tol_objective = 1e-9,
  mvnorm_algorithm = mvtnorm::Miwa(
    steps = 4097,
    checkCorr = FALSE,
    maxval = 1000),
  nloptr_opts = list(algorithm = "NLOPT_LN_SBPLX",
                     xtol_abs = 1e-10,
                     xtol_rel = 1e-9,
                     maxeval = 2000,
                     print_level = 0)
)

# This loop will create a list of parameters for the five designs compared in the
# design comparison table.
params_design_comparison <- list()
for (i in 1:2) {
  beta <- results_onestage$`$\\beta$`[i]
  # Parameters for design 2. Note that the allocation ratios come from the optimal
  # single-stage design.
  params_design_comparison[[length(params_design_comparison) + 1]] <-
    append(shared_params_twostage,
           list(
             nr = 2,
             beta = beta,
             bTP1f = -Inf, bTC1f = -Inf,
             cP1 = results_onestage$.design_object[[i]]$stagec[[1]][["P"]],
             cC1 = results_onestage$.design_object[[i]]$stagec[[1]][["C"]],
             cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
             binding_futility = FALSE
           ))
  # Parameters for design 3. No restrictions are put on the allocation ratios.
  params_design_comparison[[length(params_design_comparison) + 1]] <-
    append(shared_params_twostage,
           list(
             nr = 3,
             beta = beta,
             bTP1f = -Inf, bTC1f = -Inf,
             binding_futility = FALSE
           ))
  # Parameters for design 4. Includes non-binding futility boundaries.
  params_design_comparison[[length(params_design_comparison) + 1]] <-
    append(shared_params_twostage,
           list(
             nr = 4,
             beta = beta,
             binding_futility = FALSE
           ))
  # Parameters for design 5. Includes binding futility boundaries.
  params_design_comparison[[length(params_design_comparison) + 1]] <-
    append(shared_params_twostage,
           list(
             nr = 5,
             beta = beta,
             binding_futility = TRUE
           ))
}
# These are the design parameters for the investigation of lambda in the objective.
params_varied_lambda <- expand.grid.list(
  nr = 5,
  lambda = seq(.1, .9, .1),
  beta = c(.2, .1),
  binding_futility = TRUE) |>
  lapply(\(x)append(x, shared_params_twostage))
# These are the design parameters for the investigation of kappa in the objective.
params_varied_kappa <- expand.grid.list(
  nr = 5,
  kappa = seq(0.5, 3, by = .5),
  beta = c(.2, .1),
  binding_futility = TRUE) |>
  lapply(\(x)append(x, shared_params_twostage))
# These are the design parameters for the investigation of eta in the objective.
params_varied_eta <- expand.grid.list(
  nr = 5,
  eta = seq(0.5, 3, by = .5),
  beta = c(.2, .1),
  binding_futility = TRUE) |>
  lapply(\(x)append(x, shared_params_twostage))
# These parameters were specifically chosen, so that the optimized designs 4 and 5
# match the performance of design 2 in terms of expected sample size under the
# alternative, maximum sample size, and expected placebo group sample size, while
# achieving better performance under the null due to the futility boundaries.
params_final_comparison <- c(
  list(list(
      nr = c(4),
      lambda = .9,
      kappa = 1,
      eta = .2,
      beta = c(.2),
      binding_futility = FALSE),
    list(
      nr = c(5),
      lambda = .9,
      kappa = 1,
      eta = .2,
      beta = c(.2),
      binding_futility = TRUE)) |>
    lapply(\(x)append(x, shared_params_twostage)),
  list(list(
      nr = c(4),
      lambda = .95,
      kappa = 1.5,
      eta = .05,
      beta = c(.1),
      binding_futility = FALSE),
    list(
      nr = c(5),
      lambda = .96,
      kappa = 1.5,
      eta = .06,
      beta = c(.1),
      binding_futility = TRUE)) |>
    lapply(\(x)append(x, shared_params_twostage)))
# Combine all parameters to a single list.
all_params <- c(params_design_comparison,
                params_varied_lambda,
                params_varied_kappa,
                params_varied_eta,
                params_final_comparison)
# Helper function to perform two-step optimization. First, various starting points
# for the design parameters are tried a low-accuracy local optimization is performed
# from there. Then, the best result from this procedure is selected and a
# high-accuracy local optimization is performed.
opt_two_step <- function(...) {
  arglist <- list(...)
  orig_nloptr <- arglist$nloptr_opts
  orig_mvnorm <- arglist$mvnorm_algorithm
  arglist$nloptr_opts <- list(
    algorithm = "NLOPT_GN_MLSL_LDS",
    xtol_rel = 0.001,
    print_level = 0,
    maxeval = 10000,
    local_opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      ftol_rel = 1e-4,
      xtol_abs = 1e-3,
      xtol_rel = 1e-2,
      maxeval = 40,
      print_level = 0
    ))
  arglist$mvnorm_alogrithm <- mvtnorm::Miwa(
    steps = 128,
    checkCorr = FALSE,
    maxval = 1000)
  opt_mlsl <- do.call(optimize_design_twostage, arglist, quote = TRUE)
  arglist$nloptr_opts <- orig_nloptr
  arglist$mvnorm_alogrithm <- orig_mvnorm
  arglist$nloptr_x0 <- opt_mlsl$x_end
  opt_refined <- do.call(optimize_design_twostage, arglist, quote = TRUE)
  return(opt_refined)
}

# Set up progress bar
with_progress({
  prog <- progressor(length(all_params))
  # Optimize designs in parallel. This is where most of the calculations happen.
  results_twostage <- future_lapply(
    all_params,
    function(x) {
      prog(sprintf(
"Optimizing design %i for beta = %.2f, lambda = %.2f, kappa = %.2f and eta = %.2f.",
        x$nr,
        if (!is.null(x$beta)) x$beta else .8,
        if (!is.null(x$lambda)) x$lambda else 0,
        if (!is.null(x$kappa)) x$kappa else 0,
        if (!is.null(x$eta)) x$eta else 0
      ))
      do.call(opt_two_step, x, quote = TRUE)
    },
    future.seed = seed
  ) |>
    OptimalGoldstandardDesigns:::make_table()
})
# Combine results for single-stage designs and two-stage designs.
goldstandard_results <- rbind(results_onestage, results_twostage)
# Save the results on disk.
saveRDS(goldstandard_results, "data/goldstandard_results.rds")
# Clean up computing environment.
rm(results_onestage, results_twostage, all_params, params_design_comparison,
   params_varied_lambda, params_varied_kappa, params_varied_eta,
   params_final_comparison, shared_params_twostage, params_onestage,
   beta, i, seed, expand.grid.list, opt_two_step, prog)
