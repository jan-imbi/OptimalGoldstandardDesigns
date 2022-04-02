#' Helper function to calculate the required sample size (of the stage 1 treatment group)
#' to achieve the target power given the b2TCe 
#' 
#' This is required in designs with nonsequential futility testing, as choosing locally valid designs
#' is insufficient to gurantee type I error control.
#'
#' @template b2TCe
#' @template D
#' @include design_helper_functions.R
calc_nT1_wrt_bTC2e <- function(b2TCe, D) {
  D$b[[2]][["TC"]]["efficacy"] <- b2TCefficacy
  D$tol <- D$tol * 3 / 4
  
  if (calc_prob_reject_both(1, D) - D$type_II_error <= D$tol / 3) {
    return(1)
  } else {
    # This trick calculates (-) the sqrt of the sample size required for a power of 1-D$type_II_error to reject in the first stage.
    # This gives an upper bound on the total sample size required for a power of 1-D$type_II_error.
    sqrtn <- qmvnorm(1 - D$type_II_error,
                     mean = -as.vector(diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]) %*% c(D$b[[1]][["TP"]][["efficacy"]], D$b[[1]][["TC"]][["efficacy"]])),
                     var = diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]) %*%
                       (projection[["TP1_TC1"]] %*% D$Sigma %*% t(projection[["TP1_TC1"]])) %*% diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]),
                     tail = "upper.tail"
    )$quantile
    
    return(uniroot(function(x, D) 1 - calc_prob_reject_both(D$mu_wo_nT1[["H1"]] * sqrt(nT1), D) - D$type_II_error,
                   c(1, sqrtn^2),
                   tol = D$tol / 3,
                   extendInt = "downX",
                   D = D
    )$root)
  }
}


#' Helper function to calculate the maximal probability of rejecting the non-inferiority hypothesis
#' in the testing procedure featuring nonsequential futility, given a point hypothesis for
#' the superiority hypothesis. 
#' 
#' This is required in designs with nonsequential futility testing, as choosing locally valid designs
#' is insufficient to gurantee type I error control.
#'
#' @template b2TCe
#' @template D
calc_worst_type_I_error <- function(bTC2e, D) {
  sqrt_nT1 <- sqrt(calc_nT1_wrt_bTC2e(bTC2e, D))
  nT1_div_gamma <- sqrt_nT1 / c(D$gamma[[1]][["TP"]], D$gamma[[2]][["TP"]])
  D$b[[2]][["TC"]]["efficacy"] <- b2TCefficacy
  
  calc_alpha_wrt_muH0T <- function(muH0TP, D) {
    return(vapply(muH0TP, function(x, D){
      mu_vec <- c(x * nT1_div_gamma, 0, 0)
      calc_prob_reject_both_wrt_mu_vec(mu_vec, D)
    }, numeric(1), D = D))
  }
  
  D$tol <- D$tol * 3 / 5
  optimum <- optimize(calc_alpha_wrt_muH0T, c(0, qnorm(1 - D$tol / 3) / min(nT1_div_gamma)),
                      maximum = TRUE, D = D,
                      tol = D$tol / 3
  )
  return(optimum$objective)
}