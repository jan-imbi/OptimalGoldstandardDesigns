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

bc_d1 <- c(
  eg(print_progress=FALSE, round_n = FALSE,
     Delta = c(.2), nr = 1, beta = c(.2, .1), varP=c(1)),
  eg(print_progress=FALSE, round_n = FALSE,
     Delta = c(.4/3), nr = 1, beta = c(.2, .1), varP=c(1)),
  eg(print_progress=FALSE, round_n = FALSE,
     Delta = c(.2), nr = 1, beta = c(.2, .1), varP=c(2))
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
  varC = 1,
  round_n = FALSE,
  print_progress=FALSE,
  inner_tol_objective = 1e-05,
  mvnorm_algorithm = mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1000),
  nloptr_opts = list(algorithm = "NLOPT_LN_SBPLX", ftol_rel = 1e-04, xtol_abs = 0.001,
                     xtol_rel = 0.01, maxeval = 70, print_level = 0)
)

bc_a1<- list()
for (i in seq_along(c(.2, .1))){
  beta <- c(.2, .1)[i]
  bc_a1[[length(bc_a1) + 1]] <-
    append(
      shared_params,
      list(
        nr = 2,
        Delta = 4/3,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        cP1 = opt_onestage[[2 + i]]$stagec[[1]][["P"]], cC1 = opt_onestage[[2 + i]]$stagec[[1]][["C"]],
        cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
        binding_futility = FALSE
        )
    )
  bc_a1[[length(bc_a1) + 1]] <-
    append(
      shared_params,
      list(
        nr = 3,
        Delta = 4/3,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        binding_futility = FALSE
      )
    )
  bc_a1[[length(bc_a1) + 1]] <-
    append(
      shared_params,
      list(
        nr = 4,
        Delta = 4/3,
        varP = 1,
        beta = beta,
        binding_futility = FALSE
      )
    )
  bc_a1[[length(bc_a1) + 1]] <-
    append(
      shared_params,
      list(
        nr = 5,
        Delta = 4/3,
        varP = 1,
        beta = beta,
        binding_futility = TRUE
      )
    )
}


bc_a2<- list()
for (i in seq_along(c(.2, .1))){
  beta <- c(.2, .1)[i]
  bc_a2[[length(bc_a2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 2,
        Delta = .2,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        cP1 = opt_onestage[[4+i]]$stagec[[1]][["P"]], cC1 = opt_onestage[[4+i]]$stagec[[1]][["C"]],
        cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
        binding_futility = FALSE
        )
    )
  bc_a2[[length(bc_a2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 3,
        Delta = .2,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        binding_futility = FALSE
      )
    )
  bc_a2[[length(bc_a2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 4,
        Delta = .2,
        varP = 1,
        beta = beta,
        binding_futility = FALSE
      )
    )
  bc_a2[[length(bc_a2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 5,
        Delta = .2,
        varP = 1,
        beta = beta,
        binding_futility = TRUE
      )
    )
}

bc_atemp <- list()
for (i in seq_along(c(.2, .1))){
  beta <- c(.2, .1)[i]
  bc_atemp[[length(bc_atemp) + 1]] <-
    append(
      shared_params,
      list(
        nr = 2,
        Delta = .2,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        cP1 = opt_onestage[[0+i]]$stagec[[1]][["P"]], cC1 = opt_onestage[[0+i]]$stagec[[1]][["C"]],
        cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
        binding_futility = FALSE
        )
    )
}


bc_a3 <- lapply(eg(
  nr = 5,
  lambda = seq(1, .1, by=-.1),
  beta = c(.2, .1),
  binding_futility = FALSE
  ), function(x)append(shared_params, x))

bc_a4 <- lapply(eg(
  nr = 5,
  kappa = seq(0, 3, by=.5),
  beta = c(.2, .1),
  binding_futility = FALSE
), function(x)append(shared_params, x))

bc_a5 <- lapply(eg(
  nr = 5,
  kappa = seq(0, 3, by=1),
  lambda = seq(1, .1, by=-.1),
  beta = .2,
  binding_futility = FALSE
), function(x)append(shared_params, x))

bc_all <- c(bc_a1, bc_a2, bc_a3, bc_a4, bc_a5, bc_atemp)

opt_all <- future_lapply(
  bc_all,
  function(x)do.call(optimize_design_twostage, x, quote = TRUE),
  future.globals = future_globals,
  future.packages = future_packages,
  future.seed = NULL
)

istart <- 1
iend <- length(bc_a1)
opt_a1 <- opt_all[istart:iend]
istart <- istart + length(bc_a1)
iend <- iend + length(bc_a2)
opt_a2 <- opt_all[istart:iend]
istart <- istart + length(bc_a2)
iend <- iend + length(bc_a3)
opt_a3 <- opt_all[istart:iend]
istart <- istart + length(bc_a3)
iend <- iend + length(bc_a4)
opt_a4 <- opt_all[istart:iend]
istart <- istart + length(bc_a4)
iend <- iend + length(bc_a5)
opt_a5 <- opt_all[istart:iend]
istart <- istart + length(bc_a5)
iend <- iend + length(bc_atemp)
opt_atemp <- opt_all[istart:iend]

D_a1 <- c(
  opt_onestage[3],
  opt_a1[1:4],
  opt_onestage[4],
  opt_a1[5:8]
)
D_a2 <- c(
  opt_onestage[5],
  opt_a2[1:4],
  opt_onestage[6],
  opt_a2[5:8]
)
D_a3 <- opt_a3
D_a4 <- c(opt_atemp[2],
            opt_a4[1:7],
            opt_a1[5],
            opt_a4[8:14])
D_a5 <- c(opt_atemp[2],
          opt_a5)

saveRDS(D_a1, here("data", "D_a1.rds"))
saveRDS(D_a2, here("data", "D_a2.rds"))
saveRDS(D_a3, here("data", "D_a3.rds"))
saveRDS(D_a4, here("data", "D_a4.rds"))
saveRDS(D_a5, here("data", "D_a5.rds"))
