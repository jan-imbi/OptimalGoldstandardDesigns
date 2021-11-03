#' Optimal design parameters for design 1
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_single_stage <- function(start = NULL, D = create_Design()) {
  if (is.null(start)) {
    start <- c(1, 1)
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_single_stage,
      lb = c(1 / 300, 1 / 20),
      ub = c(3, 20),
      opts = list(
        algorithm = "NLOPT_GN_MLSL",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval = round(D$maxeval / 100)
        )
      ),
      D = D,
      return_everything = FALSE
    )

  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_single_stage,
      # eval_g_ineq = g,
      lb = c(1 / 300, 1 / 20),
      ub = c(3, 20),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}

#' Optimal design parameters for design 2
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_no_futility_fixed_c <- function(start = NULL, D = create_Design()) {
  if (is.null(start)) {
    start <- c(.5, .5)
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_no_futility_fixed_c,
      lb = c(1 / 20, 1 / 20),
      ub = c(5, 5),
      opts = list(
        algorithm = "NLOPT_GN_MLSL",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval = D$maxeval / 100
        )
      ),
      D = D,
      return_everything = FALSE
    )

  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_no_futility_fixed_c,
      lb = c(1 / 20, 1 / 20),
      ub = c(5, 5),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}

#' Optimal design parameters for design 3
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_no_futility <- function(start = NULL, D = create_Design( maxeval = 10)) {
  if (is.null(start)) {
    #         P1T C1T T2T P2T C2T
    start <- c(
      .25, 1, 1, .25, 1,
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1], # bTPE1
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1]  # b1TCE1
    )
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_no_futility,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, qnorm(1 - D$type_I_error), qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, 6, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval = D$maxeval / 100
        )
      ),
      D = D,
      return_everything = FALSE
    )
  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_no_futility,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, qnorm(1 - D$type_I_error), qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, 6, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}

#' Optimal design parameters for design 4
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_nonbinding_futility <- function(start = NULL, D = create_Design( maxeval = 10)) {
  if (is.null(start)) {
    #       P1T  C1T  T2T  P2T  C2T   bTPF1
    start <- c(
      .25, 1, 1, .25, 1, 0,
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1], # bTPE1
      0, # bTC1F
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1] # bTC1E
    )
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_nonbinding_futility,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval=D$maxeval/100
        )
      ),
      D = D,
      return_everything = FALSE
    )
  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_nonbinding_futility,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}

#' Optimal design parameters for design 6
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_closed_testing <- function(start = NULL, D = create_Design(maxeval = 10)) {
  if (is.null(start)) {
    #        P1T  C1T  T2T  P2T  C2T   bTPF1
    start <- c(
      .25, 1, 1, .25, 1, 0,
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1], # bTPE1
      0, # bTC1F
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1]
    ) # bTC1E
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_closed_testing,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval=D$maxeval / 100
        )
      ),
      D = D,
      return_everything = FALSE
    )
  if (any(opt$solution < c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error))))
    print(paste0("solution: ", opt$solution))
  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_closed_testing,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}


#' Optimal design parameters for design 5
#'
#' @param start
#' @param D
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
opt_objective_fully_sequential <- function(start = NULL, D = create_Design(maxeval = 10)) {
  if (is.null(start)) {
    #        P1T  C1T  T2T  P2T  C2T   bTPF1
    start <- c(
      .25, 1, 1, .25, 1, 0,
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1], # bTPE1
      0, # bTC0F
      gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = .5)$upper$bound[1]
    ) # bTC1E
  }
  opt <-
    nloptr(
      start,
      eval_f = objective_fully_sequential,
      # eval_g_ineq = g,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval,  # /10, (consider less maxeval because this is by far the slowest procedure)
        local_opts = list(
          algorithm = "NLOPT_LN_SBPLX",
          ftol_abs = .1,
          print_level = 0,
          maxeval=D$maxeval / 100  # /10
        )
      ),
      D = D,
      return_everything = FALSE
    )
  opt <-
    nloptr(
      opt$solution,
      eval_f = objective_fully_sequential,
      lb = c(1 / 300, 1 / 20, 1 / 20, 1 / 300, 1 / 20, -6, qnorm(1 - D$type_I_error), -6, qnorm(1 - D$type_I_error)),
      ub = c(3, 20, 20, 3, 20, qnorm(1 - D$type_I_error) - .1, 6, qnorm(1 - D$type_I_error) - .1, 6),
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = D[["tol"]] / length(start),
        print_level = 0,
        maxeval = D$maxeval  # /10
      ),
      D = D,
      return_everything = FALSE
    )
  return(opt)
}

#' Optimizes designs 1 through 6 in succession, using the previous optimum as the starting point for the next design.
#'
#' @param boundaryvals list of boundary conditions for the designs to be optimized
#' @param ncores number of cores to use for parallelization
#' @param designs vector of designs containing the numbers 1 through 6
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom nloptr nloptr
#' @import doRNG
#' @import doFuture
init_successive_optimum_designs <- function(boundaryvals = list(), ncores = NULL, designs = 1:6) {
  plan(tweak(multisession, workers = if(is.null(ncores)) length(boundaryvals) else ncores  ))
  registerDoFuture()
  erg <- foreach(
    l = boundaryvals, .errorhandling = "stop"
  ) %dorng% {
    b_target_type_I_error <- l[["type_I_error"]]
    b_target_type_II_error <- l[["type_II_error"]]
    b_Delta <- l[["Delta"]]
    b_alternative_TP <- l[["alternative_TP"]]
    b_lambda <- l[["lambda"]]
    b_kappa <- l[["kappa"]]
    b_sigma_T <- l[["sigma_T"]]
    b_sigma_P <- l[["sigma_P"]]
    b_sigma_C <- l[["sigma_C"]]
    b_tol <- l[["tol"]]
    b_maxpts <- l[["maxpts"]]
    b_maxeval <- l[["maxeval"]]

    D <- create_Design(
      type_I_error = b_target_type_I_error,
      type_II_error = b_target_type_II_error,
      Delta = b_Delta,
      alternative_TP = b_alternative_TP,
      lambda = b_lambda,
      kappa = b_kappa,
      sigma_T = b_sigma_T,
      sigma_P = b_sigma_P,
      sigma_C = b_sigma_C,
      tol = b_tol,
      maxpts = b_maxpts,
      maxeval = b_maxeval
    )

    if (1 %in% designs) {
      opt_single_stage <- opt_objective_single_stage(start = NULL, D = D)
      optD_single_stage <- objective_single_stage(opt_single_stage$solution, D, TRUE)
      optD_single_stage$opt <- opt_single_stage

      cc <- list()
      cc[[2]] <- list()
      cc[[1]][["T"]] <- 1
      cc[[1]][["P"]] <- opt_single_stage$solution[[1]]
      cc[[1]][["C"]] <- opt_single_stage$solution[[2]]
      cc[[2]][["T"]] <- 1
      cc[[2]][["P"]] <- opt_single_stage$solution[[1]]
      cc[[2]][["C"]] <- opt_single_stage$solution[[2]]
      D$cc <- cc
    } else {
      cc <- list()
      cc[[2]] <- list()
      cc[[1]][["T"]] <- 1
      cc[[1]][["P"]] <- .5
      cc[[1]][["C"]] <- 1
      cc[[2]][["T"]] <- 1
      cc[[2]][["P"]] <- .5
      cc[[2]][["C"]] <- 1
      D$cc <- cc
    }

    if (2 %in% designs) {
      opt_no_futility_fixed_c <- opt_objective_no_futility_fixed_c(start = NULL, D = D)
      optD_no_futility_fixed_c <- objective_no_futility_fixed_c(opt_no_futility_fixed_c$solution, D, TRUE)
      optD_no_futility_fixed_c$opt <- opt_no_futility_fixed_c
      wang_tp <- gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = opt_no_futility_fixed_c$solution[[1]])$upper$bound
      wang_tc <- gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = opt_no_futility_fixed_c$solution[[2]])$upper$bound
      start <- c(
        cc[[1]][["P"]],
        cc[[1]][["C"]],
        1,
        cc[[2]][["P"]],
        cc[[2]][["C"]],
        wang_tp[1],
        wang_tc[1]
      )
    } else {
      start <- c(
        cc[[1]][["P"]],
        cc[[1]][["C"]],
        1,
        cc[[2]][["P"]],
        cc[[2]][["C"]],
        2,
        2
      )
    }

    if (3 %in% designs) {
      opt_no_futility <- opt_objective_no_futility(start = start, D = D)
      optD_no_futility <- objective_no_futility(opt_no_futility$solution, D, TRUE)
      optD_no_futility$opt <- opt_no_futility

      start <- c(
        opt_no_futility$solution[[1]],
        opt_no_futility$solution[[2]],
        opt_no_futility$solution[[3]],
        opt_no_futility$solution[[4]],
        opt_no_futility$solution[[5]],
        0,
        opt_no_futility$solution[[6]],
        0,
        opt_no_futility$solution[[7]]
      )
    } else {
      start <- c(
        start[[1]],
        start[[2]],
        start[[3]],
        start[[4]],
        start[[5]],
        0,
        start[[6]],
        0,
        start[[7]]
      )
    }

    if (4 %in% designs) {
      opt_nonbinding_futility <- opt_objective_nonbinding_futility(start, D)
      optD_nonbinding_futility <- objective_nonbinding_futility(opt_nonbinding_futility$solution, D, TRUE)
      optD_nonbinding_futility$opt <- opt_nonbinding_futility

      start <- opt_nonbinding_futility$solution
      start[c(6, 8)] <- 0
    }


    if (6 %in% designs) {
      opt_closed_testing <- opt_objective_closed_testing(start, D)
      optD_closed_testing <- objective_closed_testing(opt_closed_testing$solution, D, TRUE)
      optD_closed_testing$opt <- opt_closed_testing
      start<- optD_closed_testing$solution
    }

    if (5 %in% designs) {
      # start[c(6,8)] <- 0
      opt_fully_sequential <- opt_objective_fully_sequential(start, D)
      optD_fully_sequential <- objective_fully_sequential(opt_fully_sequential$solution, D, TRUE)
      optD_fully_sequential$opt <- opt_fully_sequential
    }

    new_row <- tibble(
      target_type_I_error = b_target_type_I_error,
      target_type_II_error = b_target_type_II_error,
      Delta = b_Delta,
      alternative_TP = b_alternative_TP,
      lambda = b_lambda,
      kappa = b_kappa,
      sigma_T = b_sigma_T,
      sigma_P = b_sigma_P,
      sigma_C = b_sigma_C,
      optD_fully_sequential = if (5 %in% designs) list(optD_fully_sequential) else NULL,
      optD_closed_testing = if (6 %in% designs) list(optD_closed_testing) else NULL,
      optD_nonbinding_futility = if (4 %in% designs) list(optD_nonbinding_futility) else NULL,
      optD_no_futility = if (3 %in% designs) list(optD_no_futility) else NULL,
      optD_no_futility_fixed_c = if (2 %in% designs) list(optD_no_futility_fixed_c) else NULL,
      optD_single_stage = if (1 %in% designs) list(optD_single_stage) else NULL
    )
    new_row
  }
  return(bind_rows(erg))
}
