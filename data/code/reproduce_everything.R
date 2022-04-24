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

bc_t2<- list()
for (i in seq_along(c(.2, .1))){
  beta <- c(.2, .1)[i]
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 2,
        varP = 1,
        beta = beta,
        bTP1f = -Inf, bTC1f = -Inf,
        cP1 = opt_onestage[[i]]$stagec[[1]][["P"]], cC1 = opt_onestage[[i]]$stagec[[1]][["C"]],
        cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
        binding_futility = FALSE
      )
    )
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 3,
        varP = 1,
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
        varP = 1,
        beta = beta,
        binding_futility = FALSE
      )
    )
  bc_t2[[length(bc_t2) + 1]] <-
    append(
      shared_params,
      list(
        nr = 5,
        varP = 1,
        beta = beta,
        binding_futility = TRUE
      )
    )
}
bc_t3 <- lapply(eg(
  nr = 5,
  varP = 1,
  lambda = seq(1, .1, by=-.1),
  beta = c(.2, .1)
), function(x)append(shared_params, x))
bc_t4 <- lapply(eg(
  nr = 5,
  varP = 1,
  kappa = seq(0, 3, by=.5),
  beta = c(.2, .1)
), function(x)append(shared_params, x))
bc_t5 <- lapply(eg(
  nr = 5,
  varP = 1,
  kappa = seq(0, 3, by=1),
  lambda = seq(1, .1, by=-.1),
  beta = .2
), function(x)append(shared_params, x))

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
bc_a3 <- lapply(eg(
  nr = 5,
  varP = 1,
  lambda = seq(1, .1, by=-.1),
  beta = c(.2, .1),
  binding_futility = FALSE
), function(x)append(shared_params, x))
bc_a4 <- lapply(eg(
  nr = 5,
  varP = 1,
  kappa = seq(0, 3, by=.5),
  beta = c(.2, .1),
  binding_futility = FALSE
), function(x)append(shared_params, x))
bc_a5 <- lapply(eg(
  nr = 5,
  varP = 1,
  kappa = seq(0, 3, by=1),
  lambda = seq(1, .1, by=-.1),
  beta = .2,
  binding_futility = FALSE
), function(x)append(shared_params, x))
bc_ex4 <- append(
  shared_params,
  list(
    nr = 44,
    bTP1f = -Inf,
    bTC1f = -Inf,
    cT2 = 1,
    cP2 = quote(cP1),
    cC2 = quote(cC1),
    binding_futility = FALSE,
    beta = .1
  )
)

bc_all <- c(bc_t2, bc_t3, bc_t4, bc_t5, bc_a1, bc_a2, bc_a3, bc_a4, bc_a5, bc_ex4)
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
istart <- istart + length(opt_t5)
iend <- iend + length(opt_a1)
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
iend <- iend + length(bc_ex4)
opt_ex4 <- opt_all[istart:iend]

D_tab2 <- c(opt_onestage[1],
            opt_t2[1:4],
            opt_onestage[2],
            opt_t2[5:8])
D_tab3 <- opt_t3
D_tab4 <- c(opt_t2[2],
            opt_t4[1:7],
            opt_t2[5],
            opt_t4[8:14])
D_tab5 <- c(opt_t2[2],
            opt_t5)
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
D_ex1_tmp <- list(
  b = list(list("TP" = list("efficacy" = qnorm(1-0.025)), "TC" = list("efficacy" = qnorm(1-0.025)))),
  type_I_error = 0.025,
  type_II_error = 0.1,
  mu = list("H0" = list("TP" = 0, "TC" = 0), "H1" = list("TP" = .4, "TC" = 0 + .2)),
  var = list("T" = 1, "P" = 1, "C" = 1),
  round_n = FALSE,
  kappa = 0,
  objective = quote(sum(unlist(n)) + kappa * n[[1]][["P"]]),
  tol = 1e-7,
  mvnorm_algorithm = mvtnorm::Miwa(steps = 4097, checkCorr = FALSE, maxval = 1000),
  return_everything = TRUE,
  stagec = list(list("T" = 1, "P" = 1, "C" = 1))
)
D_ex1 <- objective_onestage(D_ex1_tmp)
D_ex2 <- opt_onestage[[2]]
D_ex3 <- opt_t2[[1]]
D_ex4 <- opt_ex4[[1]]

saveRDS(D_tab2, here("data", "D_tab2.rds"))
saveRDS(D_tab3, here("data", "D_tab3.rds"))
saveRDS(D_tab4, here("data", "D_tab4.rds"))
saveRDS(D_tab5, here("data", "D_tab5.rds"))

saveRDS(D_a1, here("data", "D_a1.rds"))
saveRDS(D_a2, here("data", "D_a2.rds"))
saveRDS(D_a3, here("data", "D_a3.rds"))
saveRDS(D_a4, here("data", "D_a4.rds"))
saveRDS(D_a5, here("data", "D_a5.rds"))

saveRDS(D_ex1, here("data", "D_ex1.rds"))
saveRDS(D_ex2, here("data", "D_ex2.rds"))
saveRDS(D_ex3, here("data", "D_ex3.rds"))
saveRDS(D_ex4, here("data", "D_ex4.rds"))

