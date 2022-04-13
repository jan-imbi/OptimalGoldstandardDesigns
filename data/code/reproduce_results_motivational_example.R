library(here)
library(mvtnorm)
library(nloptr)

invisible(sapply(paste0(here(), "/R/", list.files(here("R"))), function(x)source(x)))
future_env <- new.env()

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
D_ex2 <- optimize_design_onestage(beta = .1, round_n = FALSE, print_progress=FALSE)
D_ex3 <- optimize_design_twostage(
  bTP1f = -Inf, bTC1f = -Inf,
  cP1 = D_ex2$stagec[[1]][["P"]], cC1 = D_ex2$stagec[[1]][["C"]],
  cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
  binding_futility = FALSE,
  beta = .1, round_n = FALSE, print_progress=FALSE)
D_ex4 <- optimize_design_twostage(
  bTP1f = -Inf, bTC1f = -Inf,
  cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
  binding_futility = FALSE,
  beta = .1, round_n = FALSE, print_progress=FALSE)

saveRDS(D_ex1, here("data", "D_ex1.rds"))
saveRDS(D_ex2, here("data", "D_ex2.rds"))
saveRDS(D_ex3, here("data", "D_ex3.rds"))
saveRDS(D_ex4, here("data", "D_ex4.rds"))
