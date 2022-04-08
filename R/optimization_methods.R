#' Title
#'
#' @param cT2
#' @param cP1
#' @param cP2
#' @param cC1
#' @param cC2
#' @param bTP1f
#' @param bTP1e
#' @param i
#' @param bTC1e
#' @param alpha
#' @param beta
#' @param alternative_TP
#' @param alternative_TC
#' @param Delta
#' @param varT
#' @param varP
#' @param varC
#' @param nonbinding_futility
#' @param always_both_futility_tests
#' @param round_n
#' @param objective
#' @param inner_tol_objective
#' @param mvnorm_algorithm
#' @param nloptr_x0
#' @param nloptr_lb
#' @param nloptr_ub
#' @param nloptr_opts
#'
#' @return
#' @export
#'
#' @examples
#' 1 + 1
#'
#' @include design_helper_functions.R
optimize_design_twostage <-
  function(cT2 = 1,
           cP1 = .25,
           cP2 = .25,
           cC1 = 1,
           cC2 = 1,
           bTP1f = 0,
           bTP1e = 2.105099,
           bTC1f = 0,
           bTC1e = 2.270933,
           alpha = .025,
           beta = .2,
           alternative_TP = .4,
           alternative_TC = 0,
           Delta = .2,
           varT = 1,
           varP = 1,
           varC = 1,
           binding_futility = FALSE,
           always_both_futility_tests = TRUE,
           round_n = TRUE,
           lambda = 1,
           kappa = 0,
           nu = 0,
           objective = quote(
             sqrt(lambda) ^ 2 * ASN[["H11"]] +
               (1 - sqrt(lambda)) * sqrt(lambda) * ASN[["H10"]] +
               (1 - sqrt(lambda)) * sqrt(lambda) * ASN[["H01"]] +
               (1 - sqrt(lambda)) ^ 2 * ASN[["H00"]] +
               kappa *
               (
                 sqrt(lambda) ^ 2 * ASNP[["H11"]] +
                   (1 - sqrt(lambda)) * sqrt(lambda) * ASNP[["H10"]] +
                   (1 - sqrt(lambda)) * sqrt(lambda) * ASNP[["H01"]] +
                   (1 - sqrt(lambda)) ^ 2 * ASNP[["H00"]] +
                   nu * cumn[[2]][["P"]]
               ) +
               nu * (cumn[[2]][["T"]] + cumn[[2]][["P"]] + cumn[[2]][["C"]])
             ),
           inner_tol_objective = 1e-5,
           mvnorm_algorithm = mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3), # maxsteps = 4097
           nloptr_x0 = NULL,
           nloptr_lb = NULL,
           nloptr_ub = NULL,
           nloptr_opts = list(
             algorithm = "NLOPT_LN_SBPLX",
             ftol_rel = 1e-4,
             xtol_abs = 1e-3,
             xtol_rel = 1e-2,
             maxeval = 1000,
             print_level = 0
           ),
           print_progress = TRUE,
           ...) {
    arguments <- c(as.list(environment()))
    quoted_arguments <- arguments[sapply(arguments, is.call)]
    quoted_arguments <- quoted_arguments[names(quoted_arguments)!= "objective"]

    determine_x0 <- missing(nloptr_x0)
    determine_lb <- missing(nloptr_lb)
    determine_ub <- missing(nloptr_ub)

    if (determine_x0){
      nloptr_x0 <- list()
    }
    if (determine_lb){
      nloptr_lb <- list()
    }
    if (determine_ub){
      nloptr_ub <- list()
    }
    quotes <- list()
    quote_idx <- 1L
    if (missing(cT2)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 1
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- 1 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- 20
      quotes[[quote_idx]] <- quote(cT2 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(cP1)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- .25
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- .25 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- .25 * 20
      quotes[[quote_idx]] <- quote(cP1 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(cP2)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- .25
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- .25 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- .25 * 20
      quotes[[quote_idx]] <- quote(cP2 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(cC1)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 1
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- 1 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- 20
      quotes[[quote_idx]] <- quote(cC1 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(cC2)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 1
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- 1 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- 20
      quotes[[quote_idx]] <- quote(cC2 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(bTP1e)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 2.3
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- qnorm(1-alpha)
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- qnorm(.Machine$double.eps, lower.tail = FALSE)
      quotes[[quote_idx]] <- quote(bTP1e <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(bTP1f)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 0
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- qnorm(.Machine$double.eps, lower.tail = TRUE)
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- qnorm(1-alpha)
      quotes[[quote_idx]] <- quote(bTP1f <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(bTC1e)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 2.3
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- qnorm(1-alpha)
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- qnorm(.Machine$double.eps, lower.tail = FALSE)
      quotes[[quote_idx]] <- quote(bTC1e <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(bTC1f)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 0
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- qnorm(.Machine$double.eps, lower.tail = TRUE)
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- qnorm(1-alpha)
      quotes[[quote_idx]] <- quote(bTC1f <- x[i])
      quote_idx <- quote_idx + 1L
    }
    D_half2<- list(
      type_I_error = alpha,
      type_II_error = beta,
      mu = list("H0" = list("TP" = 0, "TC" = 0), "H1" = list("TP" = alternative_TP, "TC" = alternative_TC + Delta)),
      var = list("T" = varT, "P" = varP, "C" = varC),
      binding_futility = binding_futility,
      always_both_futility_tests = always_both_futility_tests,
      round_n = FALSE,
      lambda = lambda,
      kappa = kappa,
      nu = nu,
      objective = objective,
      tol = inner_tol_objective,
      mvnorm_algorithm = mvnorm_algorithm,
      return_everything = FALSE,
      ...
    )

    nloptr_maxeval <- if(!is.null(nloptr_opts$maxeval)) nloptr_opts$maxeval else 100
    iter <- 0

    opt_fun <- function(x){
      for (i in seq_along(x)){
        eval(quotes[[i]])
      }
      for (i in seq_along(quoted_arguments)){
        eval(quoted_arguments[[i]])
      }
      D <-
        append(list(stagec = list(
          list("T" = 1, "P" = cP1, "C" = cC1),
          list("T" = cT2, "P" = cP2, "C" = cC2)
        ),
        b = list(
          list(
            "TP" = list("efficacy" = bTP1e, "futility" = bTP1f),
            "TC" = list("efficacy" = bTC1e, "futility" = bTC1f)
          ),
          list(
            "TP" = list("efficacy" = NA_real_),
            "TC" = list("efficacy" = NA_real_)
          )
        )),
        D_half2)

      ret <- objective_twostage(D)
      if (print_progress){
        iter <<- iter + 1
        iterstring <- padd_whitespace(paste0("iteration: ", iter, "/", nloptr_maxeval))
        exprstring <- character(length = length(x))
        for (i in seq_along(x)){
          exprstring[i] <- deparse(do.call(substitute, list(quotes[[i]], list(i=i))))
        }
        exprstring <- padd_whitespace(paste0("  ",paste0(exprstring, collapse = " "), "", collapse=""))
        xstring <- padd_whitespace(paste0("    x = c(", paste0(format(x, trim=TRUE), collapse = ", "), ")", collapse = ""))
        fstring <- padd_whitespace(paste0("     f(x) = ", format(ret)))
        cat("\r", iterstring, exprstring, xstring,  fstring)
      }
      return(ret)
    }

    nloptr_x0 <- unlist(nloptr_x0)
    nloptr_lb <- unlist(nloptr_lb)
    nloptr_ub <- unlist(nloptr_ub)
    opt <- nloptr(
      x0 = nloptr_x0,
      eval_f = opt_fun,
      lb = nloptr_lb,
      ub = nloptr_ub,
      opts = nloptr_opts)

    if(print_progress)
      cat("\n")

    x <- opt$solution
    for (i in seq_along(x)){
      eval(quotes[[i]])
    }
    for (i in seq_along(quoted_arguments)){
      eval(quoted_arguments[[i]])
    }

    D_half2$round_n <- round_n
    D_half2$return_everything <- TRUE
    D <-
      append(list(stagec = list(
        list("T" = 1, "P" = cP1, "C" = cC1),
        list("T" = cT2, "P" = cP2, "C" = cC2)
      ),
      b = list(
        list(
          "TP" = list("efficacy" = bTP1e, "futility" = bTP1f),
          "TC" = list("efficacy" = bTC1e, "futility" = bTC1f)
        ),
        list(
          "TP" = list("efficacy" = NA_real_),
          "TC" = list("efficacy" = NA_real_)
        )
      )),
      D_half2)
    opt_design <- objective_twostage(D)
    opt_design$x_end <- x

    class(opt_design) <- c("TwoStageGoldStandardDesign", class(opt_design))
    return(opt_design)
  }

#' Title
#'
#' @param cT2
#' @param cP1
#' @param cP2
#' @param cC1
#' @param cC2
#' @param bTP1f
#' @param bTP1e
#' @param i
#' @param bTC1e
#' @param alpha
#' @param beta
#' @param alternative_TP
#' @param alternative_TC
#' @param Delta
#' @param varT
#' @param varP
#' @param varC
#' @param nonbinding_futility
#' @param always_both_futility_tests
#' @param round_n
#' @param objective
#' @param inner_tol_objective
#' @param mvnorm_algorithm
#' @param nloptr_x0
#' @param nloptr_lb
#' @param nloptr_ub
#' @param nloptr_opts
#'
#' @return
#' @export
#'
#' @examples
#' 1 + 1
#'
#' @include design_helper_functions.R
optimize_design_onestage <-
  function(cP1 = .25,
           cC1 = 1,
           alpha = .025,
           beta = .2,
           alternative_TP = .4,
           alternative_TC = 0,
           Delta = .2,
           varT = 1,
           varP = 1,
           varC = 1,
           round_n = TRUE,
           kappa = 0,
           objective = quote(
             sum(unlist(n)) +
               kappa * n[[1]][["P"]]
           ),
           inner_tol_objective = 1e-7,
           mvnorm_algorithm = mvtnorm::Miwa(steps = 4097, checkCorr = FALSE, maxval = 1e3), # maxsteps = 4097
           nloptr_x0 = NULL,
           nloptr_lb = NULL,
           nloptr_ub = NULL,
           nloptr_opts = list(
             algorithm = "NLOPT_LN_SBPLX",
             ftol_rel = 1e-9,
             xtol_abs = 1e-8,
             xtol_rel = 1e-7,
             maxeval = 1000,
             print_level = 0
           ),
           print_progress = TRUE,
           ...) {
    arguments <- c(as.list(environment()))
    quoted_arguments <- arguments[sapply(arguments, is.call)]
    quoted_arguments <- quoted_arguments[names(quoted_arguments)!= "objective"]

    determine_x0 <- missing(nloptr_x0)
    determine_lb <- missing(nloptr_lb)
    determine_ub <- missing(nloptr_ub)

    if (determine_x0){
      nloptr_x0 <- list()
    }
    if (determine_lb){
      nloptr_lb <- list()
    }
    if (determine_ub){
      nloptr_ub <- list()
    }
    quotes <- list()
    quote_idx <- 1L
    if (missing(cP1)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- .25
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- .25 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- .25 * 20
      quotes[[quote_idx]] <- quote(cP1 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    if (missing(cC1)){
      if (determine_x0)
        nloptr_x0[[quote_idx]] <- 1
      if (determine_lb)
        nloptr_lb[[quote_idx]] <- 1 / 20
      if (determine_ub)
        nloptr_ub[[quote_idx]] <- 20
      quotes[[quote_idx]] <- quote(cC1 <- x[i])
      quote_idx <- quote_idx + 1L
    }
    D_half2<- list(
      b = list(list("TP" = list("efficacy" = qnorm(1-alpha)), "TC" = list("efficacy" = qnorm(1-alpha)))),
      type_I_error = alpha,
      type_II_error = beta,
      mu = list("H0" = list("TP" = 0, "TC" = 0), "H1" = list("TP" = alternative_TP, "TC" = alternative_TC + Delta)),
      var = list("T" = varT, "P" = varP, "C" = varC),
      round_n = FALSE,
      kappa = kappa,
      objective = objective,
      tol = inner_tol_objective,
      mvnorm_algorithm = mvnorm_algorithm,
      return_everything = FALSE,
      ...
    )

    nloptr_maxeval <- if(!is.null(nloptr_opts$maxeval)) nloptr_opts$maxeval else 100
    iter <- 0

    opt_fun <- function(x){
      for (i in seq_along(x)){
        eval(quotes[[i]])
      }
      for (i in seq_along(quoted_arguments)){
        eval(quoted_arguments[[i]])
      }
      D <-
        append(list(stagec = list(
          list("T" = 1, "P" = cP1, "C" = cC1)
        )),
        D_half2)


      ret <- objective_onestage(D)
      if (print_progress){
        iter <<- iter + 1
        iterstring <- padd_whitespace(paste0("iteration: ", iter, "/", nloptr_maxeval))
        exprstring <- character(length = length(x))
        for (i in seq_along(x)){
          exprstring[i] <- deparse(do.call(substitute, list(quotes[[i]], list(i=i))))
        }
        exprstring <- padd_whitespace(paste0("  ",paste0(exprstring, collapse = " "), "", collapse=""))
        xstring <- padd_whitespace(paste0("    x = c(", paste0(format(x, trim = TRUE), collapse = ", "), ")", collapse = ""))
        fstring <- padd_whitespace(paste0("     f(x) = ", format(ret)))
        cat("\r", iterstring, exprstring, xstring,  fstring)
      }
      return(ret)
    }

    nloptr_x0 <- unlist(nloptr_x0)
    nloptr_lb <- unlist(nloptr_lb)
    nloptr_ub <- unlist(nloptr_ub)
    opt <- nloptr(
      x0 = nloptr_x0,
      eval_f = opt_fun,
      lb = nloptr_lb,
      ub = nloptr_ub,
      opts = nloptr_opts)

    if(print_progress)
      cat("\n")

    x <- opt$solution
    for (i in seq_along(x)){
      eval(quotes[[i]])
    }
    for (i in seq_along(quoted_arguments)){
      eval(quoted_arguments[[i]])
    }

    D_half2$round_n <- round_n
    D_half2$return_everything <- TRUE
    D <-
      append(list(stagec = list(
        list("T" = 1, "P" = cP1, "C" = cC1)
      )),
      D_half2)
    opt_design <- objective_onestage(D)
    opt_design$x_end <- x

    class(opt_design) <- c("OneStageGoldStandardDesign", class(opt_design))
    return(opt_design)
  }


optimize_with_increasing_accuracy <- function(steps = c(50, 200, 1000)){

  mvnorm_algorithm <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1000)
  nloptr_opts <- list(algorithm = "NLOPT_LN_SBPLX",
                      xtol_rel = 1e-2,
                      xtol_abs = 1e-3,
                      maxeval = 500,
                      print_level=3)

  nloptr_opts <-  list(
    algorithm = "NLOPT_GN_MLSL",
    xtol_rel = 1e-2,
    xtol_abs = 1e-3,
    print_level = 3,
    maxeval = 1000,
    local_opts = list(
      algorithm = "NLOPT_LN_SBPLX",
      ftol_rel = .1,
      print_level = 2,
      maxeval = 1000
    ))

  optimize_design_twostage(alternative_TP = .6,
                           Delta = .3,
                           round_n = FALSE,
                           inner_tol_objective = 1e-4,
                           mvnorm_algorithm = mvnorm_algorithm,
                           nloptr_opts = nloptr_opts)
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
