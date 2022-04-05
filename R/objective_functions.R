#' Objective function for two-stage designs
#'
#' @template D
#' @include design_helper_functions_design5.R
objective_twostage <- function(D)
{
  # Quote this because might need reevaluation when D$round_n == TRUE
  get_boundaries <- quote({
    
  # Calculate covariance matrix Sigma
  D$cumc <- calc_cumc(D)
  D$gamma <- calc_gamma(D)
  D$Sigma <- calc_Sigma(D)
  D$mu_wo_nT1 <- calc_mu_wo_nT1(D)
  
  # calculate b2TPefficacy boundary using the other predefined parameters
  # (This is the first implicit parameter)
  TP_local <- calc_local_rejection_boundaries("TP", D)
  alpha_TP <- TP_local$alpha
  D$b[[2]][["TP"]][["efficacy"]] <- TP_local$root
  
  # calculate bTC2e boundary using the other predefined parameters
  # (This is the second implicit parameter)
  TC_local_lower <- calc_local_rejection_boundaries("TC", D)
  
  if (D$always_both_futility_tests){
    alpha_TC <- TC_local_lower$alpha
    D$b[[2]][["TC"]][["efficacy"]] <- TC_local_lower$root
  } else{
    b2TC_lower_bound <- TC_local_lower$root
    max_ss_lb <- sum(unlist(
      calc_n_from_c(
        calc_nT1_wrt_bTC2e(b2TC_lower_bound, D), D
      )
    ))
    
    D_upper_bound <- D
    D_upper_bound$b[[1]][["TC"]][["futility"]] <- -Inf
    TC_local_upper <- calc_local_rejection_boundaries("TC", D_upper_bound)
    b2TC_upper_bound <- TC_local_upper$root
    max_ss_ub <- sum(unlist(
      calc_n_from_c(
        calc_nT1_wrt_bTC2e(b2TC_upper_bound, D), D
      )
    ))
    
    b2TC_boundary_description <- ""
    if (max_ss_ub - max_ss_lb < 1) {
      best_boundary <- list(root = b2TC_upper_bound, f.root = calc_worst_type_I_error(b2TC_upper_bound, D) - D$type_I_error)
      b2TC_boundary_description <- "boundary reduction not worth it"
    } else
      if ((lb_val <- calc_worst_type_I_error(b2TC_lower_bound, D)) - D$type_I_error <= D$tol / 5) {
        b2TC_boundary_description <- "lower boundary used"
        best_boundary <- list(root = b2TC_lower_bound, f.root = lb_val - D$type_I_error)
      } else if ((ub_val <- calc_worst_type_I_error(b2TC_upper_bound, D)) - D$type_I_error >= -D$tol / 5) {
        b2TC_boundary_description <- "upper boundary used"
        best_boundary <- list(root = b2TC_upper_bound, f.root = ub_val - D$type_I_error)
      } else {
        b2TC_boundary_description <- "optimal boundary used"
        best_boundary <- uniroot(function(x) calc_worst_type_I_error(x, D) - D$type_I_error,
                                 c(b2TC_lower_bound, b2TC_upper_bound),
                                 f.lower = lb_val, f.upper = ub_val,
                                 extendInt = "yes", tol = D$tol / 5
        )
      }
    
    D$b2TC_boundary_description <- b2TC_boundary_description
    D$b[[2]][["TC"]][["efficacy"]] <- best_boundary$root
    alpha_TC <- sum(best_boundary$f.root + D$type_I_error)
  }
  })
  eval(get_boundaries)
  
  # Calculate Design parameters, now that the implicit parameters are fixed
  nT1 <- calc_nT1_wrt_bTC2e(D$b[[2]][["TC"]][["efficacy"]], D)
  
  D$n <- n <- calc_n_from_c(nT1, D)
  if (D$round_n){
    D$n[[1]] <- lapply(D$n[[1]], ceiling)
    D$n[[2]] <- lapply(D$n[[2]], ceiling)
    D$c <- calc_c(D)
    eval(get_boundaries)
  }
  
  D$cumn <- cumn <- calc_cumn(D)
  D$mu_vec <- calc_mu_vec(D)
  
  final_state_probs <- list()
  for (hyp in c("H00", "H11", "H10", "H01")) {
    final_state_probs[[hyp]] <- calc_final_state_probs(hyp, D)
  }
  final_state_probs[["H0"]] <- final_state_probs[["H00"]]
  final_state_probs[["H1"]] <- final_state_probs[["H1"]]
  D$final_state_probs <- final_state_probs
  D$ASN <- ASN <- calc_ASN(D)
  D$ASNP <- ASNP <- calc_ASNP(D)
  lambda <- D$lambda
  kappa <- D$kappa
  nu <- D$nu
  
  objective_val <- eval(D$objective)
  
  if (D$return_everything) {
    return(D)
  } else {
    return(objective_val)
  }
}


#' Objective function for single-stage designs
#'
#' @template D
objective_singlestage <-
  function(D) {
    D$always_both_futility_tests <- FALSE
    
    # Add input parameters to Design object
    D$stagec <- list()
    D$stagec[[1]] <- list()
    D$stagec[[1]][["T"]] <- 1
    D$stagec[[1]][["P"]] <- x[1]
    D$stagec[[1]][["C"]] <- x[2]
    
    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["efficacy"]] <- qnorm(1 - D$type_I_error)
    D$b[[1]][["TC"]][["efficacy"]] <- qnorm(1 - D$type_I_error)
    
    
    # Calculate covariance matrix Sigma
    gamma <- list()
    gamma[[1]] <- list()
    for (g in c("P", "C")) {
      for (s in 1:(length(D$stagec))) {
        gamma[[s]][[paste0("T", g)]] <- sqrt(D$var[["T"]] / D$stagec[[s]][["T"]] + D$var[[g]] / D$stagec[[s]][[g]])
      }
    }
    D$gamma <- gamma
    
    Sigma <- matrix(0, ncol = 2, nrow = 2)
    Sigma[1, ] <-
      c(
        1,
        D$var[["T"]] / (D$stagec[[1]][["T"]] * gamma[[1]][["TP"]] * gamma[[1]][["TC"]])
      )
    Sigma[2, 2] <- 1
    Sigma[lower.tri(Sigma, diag = T)] <-
      t(Sigma)[lower.tri(Sigma, diag = T)]
    D$Sigma <- Sigma
    
    beta <- function(nT, D) {
      sqrt_nT <- sqrt(nT)
      mu_vec <- list()
      mu_vec[["H0"]] <-
        sqrt_nT * c(
          D$mu[["H0"]][["TP"]] / D$gamma[[1]][["TP"]],
          D$mu[["H0"]][["TC"]] / D$gamma[[1]][["TC"]]
        )
      mu_vec[["H1"]] <-
        sqrt_nT * c(
          D$mu[["H1"]][["TP"]] / D$gamma[[1]][["TP"]],
          D$mu[["H1"]][["TC"]] / D$gamma[[1]][["TC"]]
        )
      
      mu_ <- mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[2]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[2]))
      
      1 - pmvnorm_(
        mean = as.vector(mu_vec[["H1"]]),
        sigma = D$Sigma,
        lower = c(D$b[[1]][["TP"]][["efficacy"]], D$b[[1]][["TC"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
        algorithm = D$mvnorm_algorithm
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]
    }
    
    if (beta(1, D) - D$type_II_error <= D$tol / 2) {
      beta_sol <- (list(root = 1, f.root = beta(1, D) - D$type_II_error))
    } else {
      beta_sol <- (uniroot(function(x, D) beta(x, D) - D$type_II_error,
                           c(1, 10000),
                           tol = D$tol / 2,
                           extendInt = "downX",
                           D = D
      ))
    }
    
    power <- 1 - beta_sol$f.root
    nT <- beta_sol$root
    
    D$n <- list()
    D$n[[1]] <- list()
    D$n[[1]][["T"]] <- nT
    D$n[[1]][["P"]] <- nT * D$stagec[[1]][["P"]]
    D$n[[1]][["C"]] <- nT * D$stagec[[1]][["C"]]
    
    sqrt_nT <- sqrt(nT)
    mu_vec <- list()
    mu_vec[["H0"]] <-
      sqrt_nT * c(
        D$mu[["H0"]][["TP"]] / D$gamma[[1]][["TP"]],
        D$mu[["H0"]][["TC"]] / D$gamma[[1]][["TC"]]
      )
    mu_vec[["H1"]] <-
      sqrt_nT * c(
        D$mu[["H1"]][["TP"]] / D$gamma[[1]][["TP"]],
        D$mu[["H1"]][["TC"]] / D$gamma[[1]][["TC"]]
      )
    
    D$mu_vec <- mu_vec
    
    final_state_probs <- list()
    ASN <- list()
    ASNP <- list()
    
    for (hyp in c("H0", "H1")) {
      mu_ <- D$mu_vec[[hyp]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[2]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[2]))
      
      P <- list()
      P[["TP1E_TC1E"]] <- pmvnorm_(
        mean = as.vector(mu_vec[[hyp]]),
        sigma = D$Sigma,
        lower = c(D$b[[1]][["TP"]][["efficacy"]], D$b[[1]][["TC"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
        algorithm = D$mvnorm_algorithm
        # algorithm = GenzBretz(maxpts = D$mvtnorm_alg_maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]
      
      P[["TP1F"]] <- pmvnorm_(
        mean = as.vector(mu_vec[[hyp]])[1],
        sigma = D$Sigma[1, 1],
        lower = c(nInf[[1]][["TP"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]]),
        algorithm = D$mvnorm_algorithm
        # algorithm = GenzBretz(maxpts = D$mvtnorm_alg_maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]
      
      P[["TP1E_TC1F"]] <- pmvnorm_(
        mean = as.vector(mu_vec[[hyp]]),
        sigma = D$Sigma,
        lower = c(D$b[[1]][["TP"]][["efficacy"]], nInf[[1]][["TC"]]),
        upper = c(pInf[[1]][["TP"]], D$b[[1]][["TC"]][["efficacy"]]),
        algorithm = D$mvnorm_algorithm
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]
      final_state_probs[[hyp]] <- P
      ASN[[hyp]] <- sum(unlist(D$n))
      ASNP[[hyp]] <- D$n[[1]][["P"]]
    }
    
    objective_val <- eval(D$objective)
    
    if (return_everything) {
      return(D)
    } else {
      return(objective_val)
    }
  }
