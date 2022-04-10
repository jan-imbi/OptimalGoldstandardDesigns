library(future.apply)
library(here)
eg <- function(...) apply(expand.grid(...), 1, as.list)
set.seed(19102021) # Everything should deterministic, but just to be sure...

plan(list(
  tweak(multisession, workers = availableCores() - 1L) # consider using a sequential plan for low-core machines.
))

invisible(sapply(paste0(here(), "/R/", list.files(here("R"))), function(x)source(x)))
future_env <- new.env()
invisible(sapply(paste0(here(), "/R/", list.files(here("R"))), function(x)source(x, local = future_env)))
future_globals <- as.list(future_env)
future_packages <- c("mvtnorm", "nloptr")

# var_all <- 1
# alpha_all <- .025
# maxpts_all <- 4097    # consider using smaller maxpts (used for calculating mvt normal intergrals) for testing
# tol_all <- 5e-10      # consider using lower tolerance for testing
# maxeval_all <- 20000   # consider using smaller maxeval (used in the optimization routine) for testing

bc_d1 <- list(
  append(list(print_progress=FALSE, round_n = FALSE), list(nr = 1, beta = .2)),
  append(list(print_progress=FALSE, round_n = FALSE), list(nr = 1, beta = .1))
  )
opt_onestage <- 
  future_lapply(bc_d1,
                function(x)do.call(optimize_design_onestage, x, quote = TRUE),
                future.globals = future_globals,
                future.packages = future_packages,
                future.seed = NULL)

shared_params <- list(
  alpha = 0.025,
  varT = 1,
  varP = 1,
  varC = 1,
  round_n = FALSE,
  print_progress=FALSE,
  inner_tol_objective = 1e-05,
  mvnorm_algorithm = mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1000),
  nloptr_opts = list(algorithm = "NLOPT_LN_SBPLX", ftol_rel = 1e-04, xtol_abs = 0.001,
                     xtol_rel = 0.01, maxeval = 70, print_level = 0)
)

bc_t2<- list()
for (i in seq_along(c(.2, .1))){
  beta <- c(.2, .1)[i]
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 2,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        cP1 = eval(opt_onestage[[i]]$stagec[[1]][["P"]]), cC1 = eval(opt_onestage[[i]]$stagec[[1]][["C"]]),
        cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
        binding_futility = FALSE
        )
    )
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 3,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        binding_futility = FALSE
      )
    )
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 4,
        beta = beta,
        binding_futility = FALSE
      )
    )
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 5,
        beta = beta,
        binding_futility = TRUE
      )
    )
}

bc_t3 <- lapply(eg(
  nr = 5,
  lambda = seq(1, .1, by=-.1),
  beta = c(.2, .1)
  ), function(x)append(shared_params, x))

bc_t4 <- lapply(eg(
  nr = 5,
  kappa = seq(0, 3, by=.5),
  beta = c(.2, .1)
), function(x)append(shared_params, x))

bc_t5 <- lapply(eg(
  nr = 5,
  kappa = seq(0, 3, by=1),
  lambda = seq(1, .1, by=-.1),
  beta = .2
), function(x)append(shared_params, x))

bc_all <- c(bc_t2, bc_t3, bc_t4, bc_t5)

opt_all <- future_lapply(
  bc_all,
  function(x)do.call(optimize_design_twostage, x, quote = TRUE),
  future.globals = future_globals,
  future.packages = future_packages,
  future.seed = NULL
)

istart <- 1
iend <- length(bc_t2)
opt_t2 <- opt_all[istart:iend]
istart <- istart + length(bc_t2)
iend <- iend + length(bc_t3)
opt_t3 <- opt_all[istart:iend]
istart <- istart + length(bc_t3)
iend <- iend + length(bc_t4)
opt_t4 <- opt_all[istart:iend]
istart <- istart + length(bc_t4)
iend <- iend + length(bc_t5)
opt_t5 <- opt_all[istart:iend]

D_tab2 <- append(append(append(list(opt_onestage[[1]]),
                             opt_t2[1:4]),
                      list(opt_onestage[[2]])),
               opt_t2[5:8])
D_tab3 <- opt_t2
D_tab4 <- append(append(append(list(opt_t2[[2]]), opt_t4[1:7]),
                      list(opt_t2[[5]])),
               opt_t4[8:14])
D_tab5 <- append(list(opt_t2[[2]]), opt_t5)

saveRDS(value(D_tab2), here("data", "D_tab2.rds"))
saveRDS(value(D_tab3), here("data", "D_tab3.rds"))
saveRDS(value(D_tab4), here("data", "D_tab4.rds"))
saveRDS(value(D_tab5), here("data", "D_tab5.rds"))





