# library(future.apply)
# library(here)
# eg <- function(...) apply(expand.grid(...), 1, as.list)
# set.seed(19102021) # Everything should deterministic, but just to be sure...
#
# plan(list(
#   tweak(multisession, workers = availableCores() - 1L) # consider using a sequential plan for low-core machines.
# ))
#
# invisible(sapply(paste0(here(), "/R/", list.files(here("R"))), function(x)source(x)))
# future_env <- new.env()
# invisible(sapply(paste0(here(), "/R/", list.files(here("R"))), function(x)source(x, local = future_env)))
# future_globals <- as.list(future_env)
# future_packages <- c("mvtnorm", "nloptr")

D_ex1 <- list(
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
  return_everything = FALSE,
  stagec = list(list("T" = 1, "P" = 1, "C" = 1))
)
nmax_ex1 <- objective_onestage(D_ex1)
D_ex2 <- optimize_design_onestage(beta = .1, round_n = FALSE)
nT_ex2 <- D_ex2$n[[1]][["T"]] / sum(unlist(D_ex2$n[[1]]))
nP_ex2 <- D_ex2$n[[1]][["P"]] / sum(unlist(D_ex2$n[[1]]))
nC_ex2 <- D_ex2$n[[1]][["C"]] / sum(unlist(D_ex2$n[[1]]))
obj_ect2 <- D_ex2$objective_val





# saveRDS(value(D_tab2), here("data", "D_tab2.rds"))
# saveRDS(value(D_tab3), here("data", "D_tab3.rds"))
# saveRDS(value(D_tab4), here("data", "D_tab4.rds"))
# saveRDS(value(D_tab5), here("data", "D_tab5.rds"))





