#' Helper function to calculate cumulative sample size after second stage from design object
#'
#' @param D
#'
#' @return
#'
#' @examples
calc_cumn <- function(D) {
  n <- D$n
  cumn <- D$n
  cumn[[2]][["T"]] <- n[[1]][["T"]] + n[[2]][["T"]]
  cumn[[2]][["P"]] <- n[[1]][["P"]] + n[[2]][["P"]]
  cumn[[2]][["C"]] <- n[[1]][["C"]] + n[[2]][["C"]]
  return(cumn)
}

#' Helper function to calculate "cumulative allocation ratio" from design object
#'
#' @param D
#'
#' @return
#'
#' @examples
calc_cumc <- function(D) {
  cc <- D$cc
  ccc <- D$cc
  ccc[[2]][["T"]] <- cc[[1]][["T"]] + cc[[2]][["T"]]
  ccc[[2]][["P"]] <- cc[[1]][["P"]] + cc[[2]][["P"]]
  ccc[[2]][["C"]] <- cc[[1]][["C"]] + cc[[2]][["C"]]
  return(ccc)
}

#' Helper function to calculate other n's given n_{1,T} and allocation ratios from design object
#'
#' @param nT1
#' @param D
#'
#' @return
#'
#' @examples
calc_n_from_c <- function(nT1, D) {
  cc <- D$cc
  n <- list()
  n[[2]] <- list()
  n[[1]][["T"]] <- nT1
  n[[1]][["P"]] <- cc[[1]][["P"]] * n[[1]][["T"]]
  n[[1]][["C"]] <- cc[[1]][["C"]] * n[[1]][["T"]]
  n[[2]][["T"]] <- cc[[2]][["T"]] * n[[1]][["T"]]
  n[[2]][["P"]] <- cc[[2]][["P"]] * n[[1]][["T"]]
  n[[2]][["C"]] <- cc[[2]][["C"]] * n[[1]][["T"]]
  return(n)
}

#' Helper function to calculate allocation ratios from sample sizes
#'
#' @param nT1
#' @param D
#'
#' @return
#'
#' @examples
calc_c <- function(nT1, D) {
  n <- D$n
  cc <- list()
  cc[[2]] <- list()
  cc[[1]][["T"]] <- 1
  cc[[1]][["P"]] <- n[[1]][["P"]] / n[[1]][["T"]]
  cc[[1]][["C"]] <- n[[1]][["C"]] / n[[1]][["T"]]
  cc[[2]][["T"]] <- n[[2]][["T"]] / n[[1]][["T"]]
  cc[[2]][["P"]] <- n[[2]][["P"]] / n[[1]][["T"]]
  cc[[2]][["C"]] <- n[[2]][["C"]] / n[[1]][["T"]]
  return(cc)
}

#' Helper function to calculate rho from design object
#'
#' @param D
#'
#' @return
#'
#' @examples
calc_rho <- function(D) {
  sigma <- D$sigma
  ccc <- D$ccc
  rho <- list()
  rho[[1]] <- list()
  rho[[2]] <- list()
  for (g in c("P", "C")) {
    for (s in 1:(length(ccc))) {
      rho[[s]][[paste0("T", g)]] <- sqrt(sigma[["T"]] / ccc[[s]][["T"]] + sigma[[g]] / ccc[[s]][[g]])
    }
  }
  return(rho)
}

#' Helper function to covariance matrix Sigma from design object
#'
#' @param D
#'
#' @return
#' @export
#'
#' @examples
calc_Sigma <- function(D) {
  sigma <- D$sigma
  ccc <- D$ccc
  rho <- D$rho

  Sigma <- matrix(0, ncol = 4, nrow = 4)
  Sigma[1, ] <-
    c(
      1,
      rho[[2]][["TP"]] / rho[[1]][["TP"]],
      sigma[["T"]] / (ccc[[1]][["T"]] * rho[[1]][["TP"]] * rho[[1]][["TC"]]),
      sigma[["T"]] / (ccc[[2]][["T"]] * rho[[1]][["TP"]] * rho[[2]][["TC"]])
    )
  Sigma[2, c(2, 3, 4)] <-
    c(
      1,
      sigma[["T"]] / (ccc[[2]][["T"]] * rho[[2]][["TP"]] * rho[[1]][["TC"]]),
      sigma[["T"]] / (ccc[[2]][["T"]] * rho[[2]][["TP"]] * rho[[2]][["TC"]])
    )
  Sigma[3, c(3, 4)] <-
    c(
      1,
      rho[[2]][["TC"]] / rho[[1]][["TC"]]
    )
  Sigma[4, c(4)] <-
    c(1)
  Sigma[lower.tri(Sigma, diag = T)] <-
    t(Sigma)[lower.tri(Sigma, diag = T)]
  return(Sigma)
}

#' Helper function to calculate expected value for n_{1,T}=1
#'
#' This function exists for performance optimization reasons, as it saves some calculations in the long run.
#'
#' @param D
#'
#' @return
#'
#' @examples
calc_mu_wo_nT1 <- function(D) {
  rho <- D$rho
  mu <- D$mu

  mu_wo_nT1 <- list()
  mu_wo_nT1[["H0"]] <-
    c(
      mu[["H0"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H0"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H0"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H0"]][["TC"]] / rho[[2]][["TC"]]
    )
  mu_wo_nT1[["H1"]] <-
    c(
      mu[["H1"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H1"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H1"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H1"]][["TC"]] / rho[[2]][["TC"]]
    )
  return(mu_wo_nT1)
}

#' Helper function to calculate expected value from design object
#'
#' @param D
#'
#' @return
#'
#' @examples
calc_mu_vec <- function(D) {
  sqrt_nT1 <- sqrt(D$n[[1]][["T"]])
  rho <- D$rho
  mu <- D$mu

  mu_vec <- list()
  mu_vec[["H0"]] <-
    sqrt_nT1 * c(
      mu[["H0"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H0"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H0"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H0"]][["TC"]] / rho[[2]][["TC"]]
    )
  mu_vec[["H1"]] <-
    sqrt_nT1 * c(
      mu[["H1"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H1"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H1"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H1"]][["TC"]] / rho[[2]][["TC"]]
    )

  mu_vec[["H00"]] <- mu_vec[["H0"]]
  mu_vec[["H11"]] <- mu_vec[["H1"]]

  mu_vec[["H01"]] <-
    sqrt_nT1 * c(
      mu[["H0"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H0"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H1"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H1"]][["TC"]] / rho[[2]][["TC"]]
    )
  mu_vec[["H10"]] <-
    sqrt_nT1 * c(
      mu[["H1"]][["TP"]] / rho[[1]][["TP"]],
      mu[["H1"]][["TP"]] / rho[[2]][["TP"]],
      mu[["H0"]][["TC"]] / rho[[1]][["TC"]],
      mu[["H0"]][["TC"]] / rho[[2]][["TC"]]
    )

  return(mu_vec)
}

#' Helper function to create a design object
#'
#' @param type_I_error
#' @param type_II_error
#' @param Delta
#' @param alternative_TP
#' @param nonsequential_futility
#' @param lambda
#' @param kappa
#' @param sigma_T
#' @param sigma_P
#' @param sigma_C
#' @param tol
#' @param maxpts
#' @param maxeval
#'
#' @return
#' @export
#'
#' @examples
create_Design <- function(type_I_error = 0.025,
                          type_II_error = 0.1,
                          Delta = 0.3,
                          alternative_TP = 0.6,
                          lambda = 1,
                          kappa = 0,
                          sigma_T = 1,
                          sigma_P = 1,
                          sigma_C = 1,
                          tol = 1e-7,
                          maxpts = 4097,
                          maxeval = 150) {
  mu <- list()
  mu[["H0"]][["TP"]] <- 0
  mu[["H0"]][["TC"]] <- 0

  mu[["H1"]][["TP"]] <- alternative_TP
  mu[["H1"]][["TC"]] <- Delta

  sigma <- list()
  sigma[["T"]] <- sigma_T
  sigma[["P"]] <- sigma_P
  sigma[["C"]] <- sigma_C


  # Has 2 levels: [[stage]][[TP / TC]]
  A <- list()
  A[[2]] <- list()
  A[[1]][["TP"]] <- matrix(c(1, 0, 0, 0), nrow = 1, byrow = T)
  A[[2]][["TP"]] <- matrix(c(0, 1, 0, 0), nrow = 1, byrow = T)
  A[[1]][["TC"]] <- matrix(c(0, 0, 1, 0), nrow = 1, byrow = T)
  A[[2]][["TC"]] <- matrix(c(0, 0, 0, 1), nrow = 1, byrow = T)

  A_ <- list()
  A_[["TP1"]] <- A[[1]][["TP"]]
  A_[["TP12"]] <- rbind(
    A[[1]][["TP"]],
    A[[2]][["TP"]]
  )

  A_[["TC1"]] <- A[[1]][["TC"]]
  A_[["TC12"]] <- rbind(
    A[[1]][["TC"]],
    A[[2]][["TC"]]
  )
  A_[["TC2"]] <- rbind(A[[2]][["TC"]])

  A_[["TP1_TC1"]] <- rbind(
    A[[1]][["TP"]],
    A[[1]][["TC"]]
  )
  A_[["TP1_TC12"]] <- rbind(
    A[[1]][["TP"]],
    A[[1]][["TC"]],
    A[[2]][["TC"]]
  )
  A_[["TP12_TC1"]] <- rbind(
    A[[1]][["TP"]],
    A[[2]][["TP"]],
    A[[1]][["TC"]]
  )
  A_[["TP12_TC2"]] <- rbind(
    A[[1]][["TP"]],
    A[[2]][["TP"]],
    A[[2]][["TC"]]
  )

  A_[["TP12_TC12"]] <- rbind(
    A[[1]][["TP"]],
    A[[2]][["TP"]],
    A[[1]][["TC"]],
    A[[2]][["TC"]]
  )

  cc <- list()
  cc[[2]] <- list()
  cc[[1]][["T"]] <- 1
  cc[[1]][["P"]] <- .25
  cc[[1]][["C"]] <- 1
  cc[[2]][["T"]] <- 1
  cc[[2]][["P"]] <- .25
  cc[[2]][["C"]] <- 1

  return(list(
    A = A,
    A_ = A_,
    Delta = Delta,
    alternative_TP = alternative_TP,
    mu = mu,
    sigma = sigma,
    type_I_error = type_I_error,
    type_II_error = type_II_error,
    lambda = lambda,
    kappa = kappa,
    cc = cc,
    tol = tol,
    maxpts = maxpts,
    maxeval = maxeval
  ))
}


#' Function to calculate the final state probabilities
#'
#' @param hypothesis
#' @param D
#'
#' @return
#' @export
#'
#' @examples
finalStateProb <- function(hypothesis = "H0", D) {
  A_ <- D$A_
  mu_ <- D$mu_vec[[hypothesis]]
  Sigma <- D$Sigma
  b <- D$b
  nonsequential_futility <- D$nonsequential_futility

  # Hack to get MiWa algorithm to work with infinite boundaries and reasonable accuracy
  pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
  pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
  nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
  nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

  for (i in seq_len(length(b))){
    for (j in names(b[[i]])){
      for (k in names(b[[i]][[j]])){
        if (b[[i]][[j]][[k]]==Inf)
          b[[i]][[j]][[k]] <- pInf[[i]][[j]]
        if (b[[i]][[j]][[k]]==-Inf)
          b[[i]][[j]][[k]] <- nInf[[i]][[j]]
      }
    }
  }

  toladjust <- 8

  P <- list()
  P[["TP1E_TC1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
    sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / toladjust,
    #                       releps = 0)
  )[1]

  P[["TP1E_TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC12"]] %*% mu_),
    sigma = A_[["TP1_TC12"]] %*% Sigma %*% t(A_[["TP1_TC12"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["efficacy"]], pInf[[2]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / toladjust,
    #                       releps = 0)
  )[1]

  P[["TP1E_TC12F"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC12"]] %*% mu_),
    sigma = A_[["TP1_TC12"]] %*% Sigma %*% t(A_[["TP1_TC12"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], nInf[[2]][["TC"]]),
    upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / toladjust,
    #                       releps = 0)
  )[1]

  if (isTRUE(nonsequential_futility)) {
    P[["TP1F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP1"]] %*% mu_),
      sigma = A_[["TP1"]] %*% Sigma %*% t(A_[["TP1"]]),
      lower = c(nInf[[1]][["TP"]]),
      upper = c(b[[1]][["TP"]][["futility"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP1_TC1F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
      sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
      lower = c(b[[1]][["TP"]][["futility"]], nInf[[1]][["TC"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP1F_TC1F"]] <- 1 - pmvnorm_(
      mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
      sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[1]][["TC"]][["futility"]]),
      upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(
      #   maxpts = D$maxpts,
      #   abseps = D$tol / toladjust,
      #   releps = 0
      # )
    )[1]

    P[["TP12F_TC1"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC1"]] %*% mu_),
      sigma = A_[["TP12_TC1"]] %*% Sigma %*% t(A_[["TP12_TC1"]]),
      lower = c(b[[1]][["TP"]][["futility"]], nInf[[2]][["TP"]], b[[1]][["TC"]][["futility"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], b[[2]][["TP"]][["efficacy"]], pInf[[1]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP12E_TC12E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC12"]] %*% mu_),
      sigma = A_[["TP12_TC12"]] %*% Sigma %*% t(A_[["TP12_TC12"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TC"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                                abseps = D$tol / toladjust,
      #                                releps = 0)
    )[1]

    P[["TP12E_TC12F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC12"]] %*% mu_),
      sigma = A_[["TP12_TC12"]] %*% Sigma %*% t(A_[["TP12_TC12"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], nInf[[2]][["TC"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TC"]], b[[2]][["TC"]][["efficacy"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                                abseps = D$tol / toladjust,
      #                                releps = 0)
    )[1]
  } else {
    P[["TP1F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP1"]] %*% mu_),
      sigma = A_[["TP1"]] %*% Sigma %*% t(A_[["TP1"]]),
      lower = c(nInf[[1]][["TP"]]),
      upper = c(b[[1]][["TP"]][["futility"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP1E_TC1F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
      sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
      lower = c(b[[1]][["TP"]][["efficacy"]], nInf[[1]][["TC"]]),
      upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["futility"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP12F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12"]] %*% mu_),
      sigma = A_[["TP12"]] %*% Sigma %*% t(A_[["TP12"]]),
      lower = c(b[[1]][["TP"]][["futility"]], nInf[[2]][["TP"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], b[[2]][["TP"]][["efficacy"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP12E_TC2E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC2"]] %*% mu_),
      sigma = A_[["TP12_TC2"]] %*% Sigma %*% t(A_[["TP12_TC2"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]

    P[["TP12E_TC2F"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC2"]] %*% mu_),
      sigma = A_[["TP12_TC2"]] %*% Sigma %*% t(A_[["TP12_TC2"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], nInf[[2]][["TC"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], b[[2]][["TC"]][["efficacy"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / toladjust,
      #                       releps = 0)
    )[1]
  }
  return(P)
}


#' Calculate the probability to reject both hypothesis given an expected value for the test statistic and a design object
#'
#' @param mu_vec
#' @param D
#'
#' @return
#'
#' @examples
calc_prob_reject_both_wrt_mu_vec <- function(mu_vec, D) {
  A_ <- D$A_
  mu_ <- mu_vec
  Sigma <- D$Sigma
  b <- D$b
  P <- list()
  nonsequential_futility <- D$nonsequential_futility

  pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
  pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
  nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
  nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

  for (i in seq_len(length(b))){
    for (j in names(b[[i]])){
      for (k in names(b[[i]][[j]])){
        if (b[[i]][[j]][[k]]==Inf)
          b[[i]][[j]][[k]] <- pInf[[i]][[j]]
        if (b[[i]][[j]][[k]]==-Inf)
          b[[i]][[j]][[k]] <- nInf[[i]][[j]]
      }
    }
  }



  P[["TP1E_TC1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
    sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]

  P[["TP1E_TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC12"]] %*% mu_),
    sigma = A_[["TP1_TC12"]] %*% Sigma %*% t(A_[["TP1_TC12"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["efficacy"]], pInf[[2]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]

  if (isTRUE(nonsequential_futility)) {
    P[["TP12E_TC12E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC12"]] %*% mu_),
      sigma = A_[["TP12_TC12"]] %*% Sigma %*% t(A_[["TP12_TC12"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TC"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / 3,
      #                       releps = 0)
    )[1]
  } else {
    P[["TP12E_TC2E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC2"]] %*% mu_),
      sigma = A_[["TP12_TC2"]] %*% Sigma %*% t(A_[["TP12_TC2"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / 3,
      #                       releps = 0)
    )[1]
  }
  return(sum(unlist(P)))
}


#' Calculate the probability to reject both hypotheses given a choice of alternative or null point hypothesis and a design object
#'
#' @param hypothesis
#' @param D
#'
#' @return
#' @export
#'
#' @examples
calc_prob_reject_both <- function(hypothesis = "H1", D = list()) {
  A_ <- D$A_
  mu_ <- D$mu_vec[[hypothesis]]
  Sigma <- D$Sigma
  b <- D$b
  nonsequential_futility <- D$nonsequential_futility

  pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
  pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
  nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
  nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

  for (i in seq_len(length(b))){
    for (j in names(b[[i]])){
      for (k in names(b[[i]][[j]])){
        if (b[[i]][[j]][[k]]==Inf)
          b[[i]][[j]][[k]] <- pInf[[i]][[j]]
        if (b[[i]][[j]][[k]]==-Inf)
          b[[i]][[j]][[k]] <- nInf[[i]][[j]]
      }
    }
  }

  P <- list()
  P[["TP1E_TC1E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC1"]] %*% mu_),
    sigma = A_[["TP1_TC1"]] %*% Sigma %*% t(A_[["TP1_TC1"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], pInf[[1]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]

  P[["TP1E_TC12E"]] <- pmvnorm_(
    mean = as.vector(A_[["TP1_TC12"]] %*% mu_),
    sigma = A_[["TP1_TC12"]] %*% Sigma %*% t(A_[["TP1_TC12"]]),
    lower = c(b[[1]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
    upper = c(pInf[[1]][["TP"]], b[[1]][["TC"]][["efficacy"]], pInf[[2]][["TC"]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]

  if (isTRUE(nonsequential_futility)) {
    P[["TP12E_TC12E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC12"]] %*% mu_),
      sigma = A_[["TP12_TC12"]] %*% Sigma %*% t(A_[["TP12_TC12"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[1]][["TC"]][["futility"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TC"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                                abseps = D$tol / 3,
      #                                releps = 0)
    )[1]
  } else {
    P[["TP12E_TC2E"]] <- pmvnorm_(
      mean = as.vector(A_[["TP12_TC2"]] %*% mu_),
      sigma = A_[["TP12_TC2"]] %*% Sigma %*% t(A_[["TP12_TC2"]]),
      lower = c(b[[1]][["TP"]][["futility"]], b[[2]][["TP"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[2]][["TC"]]),
      algorithm = Miwa(steps = D$maxpts)
      # algorithm = GenzBretz(maxpts = D$maxpts,
      #                       abseps = D$tol / 3,
      #                       releps = 0)
    )[1]
  }
  return(sum(unlist(P)))
}

#' Calculate the average sample size
#'
#' @param P
#' @param D
#'
#' @return
#' @export
#'
#' @examples
calc_ASN <- function(P, D) {
  n <- D$n
  nonsequential_futility <- D$nonsequential_futility
  if (isTRUE(nonsequential_futility)) {
    ASN <- (P[["TP1E_TC1E"]] + P[["TP1F_TC1F"]]) * (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]]) +
      (P[["TP1E_TC12E"]] + P[["TP1E_TC12F"]]) * (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]] + n[[2]][["T"]] + n[[2]][["C"]]) +
      (P[["TP12F_TC1"]] + P[["TP12E_TC12E"]] + P[["TP12E_TC12F"]]) *
        (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]] + n[[2]][["T"]] + n[[2]][["P"]] + n[[2]][["C"]])
  } else {
    ASN <- (P[["TP1E_TC1E"]] + P[["TP1E_TC1F"]] + P[["TP1F"]]) * (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]]) +
      (P[["TP1E_TC12E"]] + P[["TP1E_TC12F"]]) * (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]] + n[[2]][["T"]] + n[[2]][["C"]]) +
      (P[["TP12F"]] + P[["TP12E_TC2E"]] + P[["TP12E_TC2F"]]) *
        (n[[1]][["T"]] + n[[1]][["P"]] + n[[1]][["C"]] + n[[2]][["T"]] + n[[2]][["P"]] + n[[2]][["C"]])
  }
  return(ASN)
}

#' Calculate the average placebo group sample size
#'
#' @param P
#' @param D
#'
#' @return
#'
#' @examples
calc_ASNP <- function(P, D) {
  n <- D$n
  nonsequential_futility <- D$nonsequential_futility

  if (isTRUE(nonsequential_futility)) {
    ASNP <- (P[["TP1E_TC1E"]] + P[["TP1F_TC1F"]]) * (n[[1]][["P"]]) +
      (P[["TP1E_TC12E"]] + P[["TP1E_TC12F"]]) * (n[[1]][["P"]]) +
      (P[["TP12F_TC1"]] + P[["TP12E_TC12E"]] + P[["TP12E_TC12F"]]) * (n[[1]][["P"]] + n[[2]][["P"]])
  } else {
    ASNP <- (P[["TP1E_TC1E"]] + P[["TP1E_TC1F"]] + P[["TP1F"]]) * (n[[1]][["P"]]) +
      (P[["TP1E_TC12E"]] + P[["TP1E_TC12F"]]) * (n[[1]][["P"]]) +
      (P[["TP12F"]] + P[["TP12E_TC2E"]] + P[["TP12E_TC2F"]]) * (n[[1]][["P"]] + n[[2]][["P"]])
  }
  return(ASNP)
}

calc_beta_wrt_nT1 <- function(nT1, D) {
  mu_vec <- D$mu_wo_nT1[["H1"]] * sqrt(nT1)
  return(1 - calc_prob_reject_both_wrt_mu_vec(mu_vec, D))
}

calc_nT1_wrt_b2TCefficacy <- function(b2TCefficacy, D) {
  D$b[[2]][["TC"]]["efficacy"] <- b2TCefficacy
  D$tol <- D$tol * 3 / 4

  if (calc_beta_wrt_nT1(1, D) - D$type_II_error <= D$tol / 3) {
    return(1)
  } else {
    # This trick calculates (-) the sqrt of the sample size required for a power of 1-D$type_II_error to reject in the first stage.
    # This gives an upper bound on the total sample size required for a power of 1-D$type_II_error.
    sqrtn <- qmvnorm(1 - D$type_II_error,
      mean = -as.vector(diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]) %*% c(D$b[[1]][["TP"]][["efficacy"]], D$b[[1]][["TC"]][["efficacy"]])),
      sigma = diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]) %*%
        (D$A_[["TP1_TC1"]] %*% D$Sigma %*% t(D$A_[["TP1_TC1"]])) %*% diag(1 / D$mu_wo_nT1[["H1"]][c(1, 3)]),
      tail = "upper.tail"
    )$quantile

    return(uniroot(function(x, D) calc_beta_wrt_nT1(x, D) - D$type_II_error,
      c(1, sqrtn^2),
      tol = D$tol / 3,
      extendInt = "downX",
      D = D
    )$root)
  }
}

calc_alpha_wrt_muH0TP <- function(muH0TP, D) {
  mu_vec <- c(D$tmp * muH0TP, 0, 0)
  prob <- calc_prob_reject_both_wrt_mu_vec(mu_vec, D)
  return(prob)
}
calc_alpha_wrt_muH0TP_vectorized <- function(muH0TP, D) {
  return(vapply(muH0TP, calc_alpha_wrt_muH0TP, numeric(1), D = D))
}

worst_type_I_error <- function(b2TCefficacy, D) {
  sqrt_nT1 <- sqrt(calc_nT1_wrt_b2TCefficacy(b2TCefficacy, D))
  D$tmp <- sqrt_nT1 / c(D$rho[[1]][["TP"]], D$rho[[2]][["TP"]])
  D$b[[2]][["TC"]]["efficacy"] <- b2TCefficacy



  D$tol <- D$tol * 3 / 5
  optimum <- optimize(calc_alpha_wrt_muH0TP_vectorized, c(0, qnorm(1 - D$tol / 3) / D$tmp),
    maximum = T, D = D,
    tol = D$tol / 3
  )
  return(optimum$objective)
}

calc_local_rejection_boundaries <- function(groups = "TP", D) {
  mu_vec <- c(0, 0, 0, 0)

  b <- D$b

  pInf <- qnorm(D$tol * 1e-2, mean = mu_vec, lower.tail = FALSE)
  pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
  nInf <- qnorm(D$tol * 1e-2, mean = mu_vec, lower.tail = TRUE)
  nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

  for (i in seq_len(length(b))){
    for (j in names(b[[i]])){
      for (k in names(b[[i]][[j]])){
        if (!is.na(b[[i]][[j]][[k]])){
          if (b[[i]][[j]][[k]]==Inf)
            b[[i]][[j]][[k]] <- pInf[[i]][[j]]
          if (b[[i]][[j]][[k]]==-Inf)
            b[[i]][[j]][[k]] <- nInf[[i]][[j]]
        }
      }
    }
  }

  P <- list()
  P[[paste0(groups, "1E")]] <- pmvnorm_(
    mean = as.vector(D$A_[[paste0(groups, "1")]] %*% mu_vec),
    sigma = D$A_[[paste0(groups, "1")]] %*% D$Sigma %*% t(D$A_[[paste0(groups, "1")]]),
    lower = b[[1]][[groups]][["efficacy"]],
    upper = pInf[[1]][[groups]],
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]
  sgn_low <- sign(pmvnorm_(
    mean = as.vector(D$A_[[paste0(groups, "12")]] %*% mu_vec),
    sigma = D$A_[[paste0(groups, "12")]] %*% D$Sigma %*% t(D$A_[[paste0(groups, "12")]]),
    lower = c(b[[1]][[groups]][["futility"]], qnorm(1 - D$type_I_error)),
    upper = c(b[[1]][[groups]][["efficacy"]], pInf[[1]][[groups]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1] + P[[paste0(groups, "1E")]] - D$type_I_error)
  sgn_high <- sign(pmvnorm_(
    mean = as.vector(D$A_[[paste0(groups, "12")]] %*% mu_vec),
    sigma = D$A_[[paste0(groups, "12")]] %*% D$Sigma %*% t(D$A_[[paste0(groups, "12")]]),
    lower = c(b[[1]][[groups]][["futility"]], qnorm(1 - D$tol / 3)),
    upper = c(b[[1]][[groups]][["efficacy"]], pInf[[1]][[groups]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1] + P[[paste0(groups, "1E")]] - D$type_I_error)

  if (sgn_high >= -D$tol / 3) {
    b2 <- list()
    b2$root <- Inf
  } else if (sgn_low <= D$tol / 3) {
    b2 <- list()
    b2$root <- qnorm(1 - D$type_I_error)
  } else {
    b2 <- uniroot(
      function(x) {
        pmvnorm_(
          mean = as.vector(D$A_[[paste0(groups, "12")]] %*% mu_vec),
          sigma = D$A_[[paste0(groups, "12")]] %*% D$Sigma %*% t(D$A_[[paste0(groups, "12")]]),
          lower = c(b[[1]][[groups]][["futility"]], x),
          upper = c(b[[1]][[groups]][["efficacy"]], pInf[[1]][[groups]]),
          algorithm = Miwa(steps = D$maxpts)
          # algorithm = GenzBretz(maxpts = D$maxpts,
          #                       abseps = D$tol / 3,
          #                       releps = 0)
        )[1] + P[[paste0(groups, "1E")]] - D$type_I_error
      },
      c(qnorm(1 - D$type_I_error), qnorm(1 - D$tol / 3)),
      tol = D$tol / 3, extendInt = "downX"
    )
  }
  b[[2]][[groups]][["efficacy"]] <- b2$root
  P[[paste0(groups, "12E")]] <- pmvnorm_(
    mean = as.vector(D$A_[[paste0(groups, "12")]] %*% mu_vec),
    sigma = D$A_[[paste0(groups, "12")]] %*% D$Sigma %*% t(D$A_[[paste0(groups, "12")]]),
    lower = c(b[[1]][[groups]][["futility"]], b[[2]][[groups]][["efficacy"]]),
    upper = c(b[[1]][[groups]][["efficacy"]], pInf[[1]][[groups]]),
    algorithm = Miwa(steps = D$maxpts)
    # algorithm = GenzBretz(maxpts = D$maxpts,
    #                       abseps = D$tol / 3,
    #                       releps = 0)
  )[1]

  return(list(
    root = b2$root,
    alpha = sum(c(P[[paste0(groups, "1E")]], P[[paste0(groups, "12E")]]))
  ))
}


#' Design 6
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_closed_testing <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- TRUE

    # Add input parameters to Design object
    D$cc <- list()
    D$cc[[2]] <- list()
    D$cc[[1]][["T"]] <- 1
    D$cc[[1]][["P"]] <- x[1]
    D$cc[[1]][["C"]] <- x[2]
    D$cc[[2]][["T"]] <- x[3]
    D$cc[[2]][["P"]] <- x[4]
    D$cc[[2]][["C"]] <- x[5]

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["futility"]] <- x[6]
    D$b[[1]][["TP"]][["efficacy"]] <- x[7]
    D$b[[1]][["TC"]][["futility"]] <- x[8]
    D$b[[1]][["TC"]][["efficacy"]] <- x[9]
    # b[[2]][["TP"]][["futility"]] <- x[13]
    D$b[[2]][["TP"]][["efficacy"]] <- NaN
    # b[[2]][["TC"]][["futility"]] <- x[14]
    D$b[[2]][["TC"]][["efficacy"]] <- NaN

    # Calculate covariance matrix Sigma
    D$ccc <- calc_cumc(D)
    D$rho <- calc_rho(D)
    D$Sigma <- calc_Sigma(D)
    D$mu_wo_nT1 <- calc_mu_wo_nT1(D)

    # calculate b2TPefficacy boundary using the other predefined parameters
    # (This is the first implicit parameter)
    TP_local <- calc_local_rejection_boundaries("TP", D)
    alpha_TP <- TP_local$alpha
    D$b[[2]][["TP"]][["efficacy"]] <- TP_local$root

    # calculate b2TCefficacy boundary using the other predefined parameters
    # (This is the second implicit parameter)
    TC_local_lower <- calc_local_rejection_boundaries("TC", D)
    alpha_TC <- TC_local_lower$alpha
    D$b[[2]][["TC"]][["efficacy"]] <- TC_local_lower$root

	# Calculate Design parameters, now that the implicit parameters are fixed



    nT1 <- calc_nT1_wrt_b2TCefficacy(D$b[[2]][["TC"]][["efficacy"]], D)
    D$n <- calc_n_from_c(nT1, D)
    D$cumn <- calc_cumn(D)
    D$mu_vec <- calc_mu_vec(D)

    finalStateProbs <- list()
    ASN <- list()
    ASNP <- list()

    for (hyp in c("H00", "H11", "H10", "H01")) {
      finalStateProbs[[hyp]] <- finalStateProb(hyp, D)
      ASN[[hyp]] <- calc_ASN(finalStateProbs[[hyp]], D)
      ASNP[[hyp]] <- calc_ASNP(finalStateProbs[[hyp]], D)
    }

    finalStateProbs[["H0"]] <- finalStateProbs[["H00"]]
    ASN[["H0"]] <- ASN[["H00"]]
    ASNP[["H0"]] <- ASNP[["H00"]]
    finalStateProbs[["H1"]] <- finalStateProbs[["H11"]]
    ASN[["H1"]] <- ASN[["H11"]]
    ASNP[["H1"]] <- ASNP[["H11"]]

    lambda_ASN <- sqrt(D$lambda)^2 * ASN[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASN[["H00"]]

    lambda_ASNP <- sqrt(D$lambda)^2 * ASNP[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASNP[["H00"]]


    # lambda_ASN <- D$lambda * ASN[["H11"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H00"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H01"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H10"]]
    #
    # lambda_ASNP <- D$lambda * ASNP[["H11"]] +
    #   (1 - D$lambda)/3 * ASNP[["H00"]] +
    #   (1 - D$lambda)/3 * ASNP[["H01"]] +
    #   (1 - D$lambda)/3 * ASNP[["H10"]]

    objective_val <- lambda_ASN + D$kappa * lambda_ASNP

    if (return_everything) {
      mu_ <- D$mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$A_[["TP1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP1"]] %*% D$Sigma %*% t(D$A_[["TP1"]]),
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1] + pmvnorm_(
        mean = as.vector(D$A_[["TP12_TC1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP12_TC1"]] %*% D$Sigma %*% t(D$A_[["TP12_TC1"]]),
        lower = c(D$b[[1]][["TP"]][["futility"]], D$b[[2]][["TP"]][["efficacy"]], D$b[[1]][["TP"]][["futility"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]

      D$power <- calc_prob_reject_both("H1", D)
      D$min_conditional_power <- calc_conditional_power(
        D$b[[1]][["TP"]][["futility"]] + 1e-10, D$b[[1]][["TC"]][["futility"]] + 1e-10,
        D$mu_vec[["H1"]], D$Sigma, D$b, D$nonsequential_futility
      )

      D$alpha_TP <- alpha_TP
      D$alpha_TC <- alpha_TC
      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }


#' Design 5
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_fully_sequential <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- FALSE

    # Add input parameters to Design object
    D$cc <- list()
    D$cc[[2]] <- list()
    D$cc[[1]][["T"]] <- 1
    D$cc[[1]][["P"]] <- x[1]
    D$cc[[1]][["C"]] <- x[2]
    D$cc[[2]][["T"]] <- x[3]
    D$cc[[2]][["P"]] <- x[4]
    D$cc[[2]][["C"]] <- x[5]

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["futility"]] <- x[6]
    D$b[[1]][["TP"]][["efficacy"]] <- x[7]
    D$b[[1]][["TC"]][["futility"]] <- x[8]
    D$b[[1]][["TC"]][["efficacy"]] <- x[9]
    # b[[2]][["TP"]][["futility"]] <- x[13]
    D$b[[2]][["TP"]][["efficacy"]] <- NaN
    # b[[2]][["TC"]][["futility"]] <- x[14]
    D$b[[2]][["TC"]][["efficacy"]] <- NaN

    # Calculate covariance matrix Sigma
    D$ccc <- calc_cumc(D)
    D$rho <- calc_rho(D)
    D$Sigma <- calc_Sigma(D)
    D$mu_wo_nT1 <- calc_mu_wo_nT1(D)

    # calculate b2TPefficacy boundary using the other predefined parameters
    # (This is the first implicit parameter)
    TP_local <- calc_local_rejection_boundaries("TP", D)
    alpha_TP <- TP_local$alpha
    D$b[[2]][["TP"]][["efficacy"]] <- TP_local$root

    # calculate b2TCefficacy boundary using the other predefined parameters
    # (This is the second implicit parameter)
    TC_local_lower <- calc_local_rejection_boundaries("TC", D)
    b2TC_lower_bound <- TC_local_lower$root
    max_ss_lb <- sum(unlist(
      calc_n_from_c(
        calc_nT1_wrt_b2TCefficacy(b2TC_lower_bound, D), D
      )
    ))

    D_upper_bound <- D
    D_upper_bound$b[[1]][["TC"]][["futility"]] <- -Inf
    TC_local_upper <- calc_local_rejection_boundaries("TC", D_upper_bound)
    b2TC_upper_bound <- TC_local_upper$root
    max_ss_ub <- sum(unlist(
      calc_n_from_c(
        calc_nT1_wrt_b2TCefficacy(b2TC_upper_bound, D), D
      )
    ))

    # print(b2TC_upper_bound - b2TC_lower_bound)
    b2TC_boundary_description <- ""

    if (max_ss_ub - max_ss_lb < 1) {
      best_boundary <- list(root = b2TC_upper_bound, f.root = worst_type_I_error(b2TC_upper_bound, D) - D$type_I_error)
      # print("boundary reduction not worth it, upper boundary used")
      b2TC_boundary_description <- "boundary reduction not worth it"
    } else
    if ((lb_val <- worst_type_I_error(b2TC_lower_bound, D)) - D$type_I_error <= D$tol / 5) {
      b2TC_boundary_description <- "lower boundary used"
      best_boundary <- list(root = b2TC_lower_bound, f.root = lb_val - D$type_I_error)
    } else if ((ub_val <- worst_type_I_error(b2TC_upper_bound, D)) - D$type_I_error >= -D$tol / 5) {
      b2TC_boundary_description <- "upper boundary used"
      best_boundary <- list(root = b2TC_upper_bound, f.root = ub_val - D$type_I_error)
    } else {
      b2TC_boundary_description <- "optimal boundary used"
      best_boundary <- uniroot(function(x) worst_type_I_error(x, D) - D$type_I_error,
        c(b2TC_lower_bound, b2TC_upper_bound),
        f.lower = lb_val, f.upper = ub_val,
        extendInt = "yes", tol = D$tol / 5
      )
    }

    D$local_design_was_used <- b2TC_boundary_description
    D$b[[2]][["TC"]][["efficacy"]] <- best_boundary$root
    alpha_TC <- sum(best_boundary$f.root + D$type_I_error)

    # Calculate Design parameters, now that the implicit parameters are fixed
    nT1 <- calc_nT1_wrt_b2TCefficacy(D$b[[2]][["TC"]][["efficacy"]], D)
    D$n <- calc_n_from_c(nT1, D)
    D$cumn <- calc_cumn(D)
    D$mu_vec <- calc_mu_vec(D)

    finalStateProbs <- list()
    ASN <- list()
    ASNP <- list()

    for (hyp in c("H00", "H11", "H10", "H01")) {
      finalStateProbs[[hyp]] <- finalStateProb(hyp, D)
      ASN[[hyp]] <- calc_ASN(finalStateProbs[[hyp]], D)
      ASNP[[hyp]] <- calc_ASNP(finalStateProbs[[hyp]], D)
    }

    finalStateProbs[["H0"]] <- finalStateProbs[["H00"]]
    ASN[["H0"]] <- ASN[["H00"]]
    ASNP[["H0"]] <- ASNP[["H00"]]
    finalStateProbs[["H1"]] <- finalStateProbs[["H11"]]
    ASN[["H1"]] <- ASN[["H11"]]
    ASNP[["H1"]] <- ASNP[["H11"]]

    lambda_ASN <- sqrt(D$lambda)^2 * ASN[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASN[["H00"]]

    lambda_ASNP <- sqrt(D$lambda)^2 * ASNP[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASNP[["H00"]]


    # lambda_ASN <- D$lambda * ASN[["H11"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H00"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H01"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H10"]]
    #
    # lambda_ASNP <- D$lambda * ASNP[["H11"]] +
    #   (1 - D$lambda)/3 * ASNP[["H00"]] +
    #   (1 - D$lambda)/3 * ASNP[["H01"]] +
    #   (1 - D$lambda)/3 * ASNP[["H10"]]

    objective_val <- lambda_ASN + D$kappa * lambda_ASNP

    if (return_everything) {
      mu_ <- D$mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$A_[["TP1_TC1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP1_TC1"]] %*% D$Sigma %*% t(D$A_[["TP1_TC1"]]),
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1] + pmvnorm_(
        mean = as.vector(D$A_[["TP12"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP12"]] %*% D$Sigma %*% t(D$A_[["TP12"]]),
        lower = c(D$b[[1]][["TP"]][["futility"]], D$b[[2]][["TP"]][["efficacy"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]

      D$power <- calc_prob_reject_both("H1", D)
      D$min_conditional_power <- calc_conditional_power(
        D$b[[1]][["TP"]][["futility"]] + 1e-10, D$b[[1]][["TC"]][["futility"]] + 1e-10,
        D$mu_vec[["H1"]], D$Sigma, D$b, D$nonsequential_futility
      )

      D$alpha_TP <- alpha_TP
      D$alpha_TC <- alpha_TC
      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }

#' Design 4
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_nonbinding_futility <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- TRUE

    # Add input parameters to Design object
    D$cc <- list()
    D$cc[[2]] <- list()
    D$cc[[1]][["T"]] <- 1
    D$cc[[1]][["P"]] <- x[1]
    D$cc[[1]][["C"]] <- x[2]
    D$cc[[2]][["T"]] <- x[3]
    D$cc[[2]][["P"]] <- x[4]
    D$cc[[2]][["C"]] <- x[5]

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["futility"]] <- -Inf
    D$b[[1]][["TP"]][["efficacy"]] <- x[7]
    D$b[[1]][["TC"]][["futility"]] <- -Inf
    D$b[[1]][["TC"]][["efficacy"]] <- x[9]
    # b[[2]][["TP"]][["futility"]] <- x[13]
    D$b[[2]][["TP"]][["efficacy"]] <- NaN
    # b[[2]][["TC"]][["futility"]] <- x[14]
    D$b[[2]][["TC"]][["efficacy"]] <- NaN

    # Calculate covariance matrix Sigma
    D$ccc <- calc_cumc(D)
    D$rho <- calc_rho(D)
    D$Sigma <- calc_Sigma(D)
    D$mu_wo_nT1 <- calc_mu_wo_nT1(D)

    # calculate b2TPefficacy boundary using the other predefined parameters
    # (This is the first implicit parameter)
    TP_local <- calc_local_rejection_boundaries("TP", D)
    alpha_TP <- TP_local$alpha
    D$b[[2]][["TP"]][["efficacy"]] <- TP_local$root


    # calculate b2TCefficacy boundary using the other predefined parameters
    # (This is the second implicit parameter)
    TC_local_lower <- calc_local_rejection_boundaries("TC", D)
    alpha_TC <- TC_local_lower$alpha
    D$b[[2]][["TC"]][["efficacy"]] <- TC_local_lower$root

    # Calculate Design parameters, now that the implicit parameters are fixed
    D$b[[1]][["TP"]][["futility"]] <- x[6]
    D$b[[1]][["TC"]][["futility"]] <- x[8]

    nT1 <- calc_nT1_wrt_b2TCefficacy(D$b[[2]][["TC"]][["efficacy"]], D)
    D$n <- calc_n_from_c(nT1, D)
    D$cumn <- calc_cumn(D)
    D$mu_vec <- calc_mu_vec(D)

    finalStateProbs <- list()
    ASN <- list()
    ASNP <- list()

    for (hyp in c("H00", "H11", "H10", "H01")) {
      finalStateProbs[[hyp]] <- finalStateProb(hyp, D)
      ASN[[hyp]] <- calc_ASN(finalStateProbs[[hyp]], D)
      ASNP[[hyp]] <- calc_ASNP(finalStateProbs[[hyp]], D)
    }

    finalStateProbs[["H0"]] <- finalStateProbs[["H00"]]
    ASN[["H0"]] <- ASN[["H00"]]
    ASNP[["H0"]] <- ASNP[["H00"]]
    finalStateProbs[["H1"]] <- finalStateProbs[["H11"]]
    ASN[["H1"]] <- ASN[["H11"]]
    ASNP[["H1"]] <- ASNP[["H11"]]

    lambda_ASN <- sqrt(D$lambda)^2 * ASN[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASN[["H00"]]

    lambda_ASNP <- sqrt(D$lambda)^2 * ASNP[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASNP[["H00"]]


    # lambda_ASN <- D$lambda * ASN[["H11"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H00"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H01"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H10"]]
    #
    # lambda_ASNP <- D$lambda * ASNP[["H11"]] +
    #   (1 - D$lambda)/3 * ASNP[["H00"]] +
    #   (1 - D$lambda)/3 * ASNP[["H01"]] +
    #   (1 - D$lambda)/3 * ASNP[["H10"]]
    objective_val <- lambda_ASN + D$kappa * lambda_ASNP

    if (return_everything) {
      mu_ <- D$mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$A_[["TP1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP1"]] %*% D$Sigma %*% t(D$A_[["TP1"]]),
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1] + pmvnorm_(
        mean = as.vector(D$A_[["TP12_TC1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP12_TC1"]] %*% D$Sigma %*% t(D$A_[["TP12_TC1"]]),
        lower = c(D$b[[1]][["TP"]][["futility"]], D$b[[2]][["TP"]][["efficacy"]], D$b[[1]][["TP"]][["futility"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]], pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]

      D$power <- calc_prob_reject_both("H1", D)
      D$min_conditional_power <- calc_conditional_power(
        D$b[[1]][["TP"]][["futility"]] + 1e-10, D$b[[1]][["TC"]][["futility"]] + 1e-10,
        D$mu_vec[["H1"]], D$Sigma, D$b, D$nonsequential_futility
      )

      D$alpha_TP <- alpha_TP
      D$alpha_TC <- alpha_TC
      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }

#' Design 3
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_no_futility <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- FALSE

    # Add input parameters to Design object
    D$cc <- list()
    D$cc[[2]] <- list()
    D$cc[[1]][["T"]] <- 1
    D$cc[[1]][["P"]] <- x[1]
    D$cc[[1]][["C"]] <- x[2]
    D$cc[[2]][["T"]] <- x[3]
    D$cc[[2]][["P"]] <- x[4]
    D$cc[[2]][["C"]] <- x[5]

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["futility"]] <- -Inf
    D$b[[1]][["TP"]][["efficacy"]] <- x[6]
    D$b[[1]][["TC"]][["futility"]] <- -Inf
    D$b[[1]][["TC"]][["efficacy"]] <- x[7]
    # b[[2]][["TP"]][["futility"]] <- x[13]
    D$b[[2]][["TP"]][["efficacy"]] <- NaN
    # b[[2]][["TC"]][["futility"]] <- x[14]
    D$b[[2]][["TC"]][["efficacy"]] <- NaN

    # Calculate covariance matrix Sigma
    D$ccc <- calc_cumc(D)
    D$rho <- calc_rho(D)
    D$Sigma <- calc_Sigma(D)
    D$mu_wo_nT1 <- calc_mu_wo_nT1(D)

    # calculate b2TPefficacy boundary using the other predefined parameters
    # (This is the first implicit parameter)
    TP_local <- calc_local_rejection_boundaries("TP", D)
    alpha_TP <- TP_local$alpha
    D$b[[2]][["TP"]][["efficacy"]] <- TP_local$root

    # calculate b2TCefficacy boundary using the other predefined parameters
    # (This is the second implicit parameter)

    TC_local_lower <- calc_local_rejection_boundaries("TC", D)
    D$b[[2]][["TC"]][["efficacy"]] <- TC_local_lower$root
    alpha_TC <- TC_local_lower$alpha


    # Calculate Design parameters, now that the implicit parameters are fixed
    nT1 <- calc_nT1_wrt_b2TCefficacy(D$b[[2]][["TC"]][["efficacy"]], D)
    D$n <- calc_n_from_c(nT1, D)
    D$cumn <- calc_cumn(D)
    D$mu_vec <- calc_mu_vec(D)

    finalStateProbs <- list()
    ASN <- list()
    ASNP <- list()

    for (hyp in c("H00", "H11", "H10", "H01")) {
      finalStateProbs[[hyp]] <- finalStateProb(hyp, D)
      ASN[[hyp]] <- calc_ASN(finalStateProbs[[hyp]], D)
      ASNP[[hyp]] <- calc_ASNP(finalStateProbs[[hyp]], D)
    }

    finalStateProbs[["H0"]] <- finalStateProbs[["H00"]]
    ASN[["H0"]] <- ASN[["H00"]]
    ASNP[["H0"]] <- ASNP[["H00"]]
    finalStateProbs[["H1"]] <- finalStateProbs[["H11"]]
    ASN[["H1"]] <- ASN[["H11"]]
    ASNP[["H1"]] <- ASNP[["H11"]]

    lambda_ASN <- sqrt(D$lambda)^2 * ASN[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASN[["H00"]]

    lambda_ASNP <- sqrt(D$lambda)^2 * ASNP[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASNP[["H00"]]


    # lambda_ASN <- D$lambda * ASN[["H11"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H00"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H01"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H10"]]
    #
    # lambda_ASNP <- D$lambda * ASNP[["H11"]] +
    #   (1 - D$lambda)/3 * ASNP[["H00"]] +
    #   (1 - D$lambda)/3 * ASNP[["H01"]] +
    #   (1 - D$lambda)/3 * ASNP[["H10"]]

    objective_val <- lambda_ASN + D$kappa * lambda_ASNP


    if (return_everything) {
      mu_ <- D$mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$A_[["TP1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP1"]] %*% D$Sigma %*% t(D$A_[["TP1"]]),
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1] + pmvnorm_(
        mean = as.vector(D$A_[["TP12"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP12"]] %*% D$Sigma %*% t(D$A_[["TP12"]]),
        lower = c(nInf[[1]][["TP"]], D$b[[2]][["TP"]][["efficacy"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]

      D$power <- calc_prob_reject_both("H1", D)
      D$min_conditional_power <- 0


      D$alpha_TP <- alpha_TP
      D$alpha_TC <- alpha_TC

      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }

#' Design 2
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_no_futility_fixed_c <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- FALSE

    wang_tp <- gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = x[1], tol = D$tol / 2)$upper$bound
    wang_tc <- gsDesign::gsDesign(2, 1, sfu = "WT", sfupar = x[2], tol = D$tol / 2)$upper$bound

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["futility"]] <- -Inf
    D$b[[1]][["TP"]][["efficacy"]] <- wang_tp[1]
    D$b[[1]][["TC"]][["futility"]] <- -Inf
    D$b[[1]][["TC"]][["efficacy"]] <- wang_tc[1]
    D$b[[2]][["TP"]][["efficacy"]] <- wang_tp[2]
    D$b[[2]][["TC"]][["efficacy"]] <- wang_tc[2]

    # Calculate covariance matrix Sigma
    D$ccc <- calc_cumc(D)
    D$rho <- calc_rho(D)
    D$Sigma <- calc_Sigma(D)
    D$mu_wo_nT1 <- calc_mu_wo_nT1(D)

    # Calculate Design parameters, now that the implicit parameters are fixed
    nT1 <- calc_nT1_wrt_b2TCefficacy(D$b[[2]][["TC"]][["efficacy"]], D)
    D$n <- calc_n_from_c(nT1, D)
    D$cumn <- calc_cumn(D)
    D$mu_vec <- calc_mu_vec(D)

    finalStateProbs <- list()
    ASN <- list()
    ASNP <- list()

    for (hyp in c("H00", "H11", "H10", "H01")) {
      finalStateProbs[[hyp]] <- finalStateProb(hyp, D)
      ASN[[hyp]] <- calc_ASN(finalStateProbs[[hyp]], D)
      ASNP[[hyp]] <- calc_ASNP(finalStateProbs[[hyp]], D)
    }

    finalStateProbs[["H0"]] <- finalStateProbs[["H00"]]
    ASN[["H0"]] <- ASN[["H00"]]
    ASNP[["H0"]] <- ASNP[["H00"]]
    finalStateProbs[["H1"]] <- finalStateProbs[["H11"]]
    ASN[["H1"]] <- ASN[["H11"]]
    ASNP[["H1"]] <- ASNP[["H11"]]

    lambda_ASN <- sqrt(D$lambda)^2 * ASN[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASN[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASN[["H00"]]

    lambda_ASNP <- sqrt(D$lambda)^2 * ASNP[["H11"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H10"]] +
      (1 - sqrt(D$lambda))*sqrt(D$lambda) * ASNP[["H01"]] +
      (1 - sqrt(D$lambda))^2 * ASNP[["H00"]]


    # lambda_ASN <- D$lambda * ASN[["H11"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H00"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H01"]] +
    #   (1 - sqrt(D$lambda))/ 3 * ASN[["H10"]]
    #
    # lambda_ASNP <- D$lambda * ASNP[["H11"]] +
    #   (1 - D$lambda)/3 * ASNP[["H00"]] +
    #   (1 - D$lambda)/3 * ASNP[["H01"]] +
    #   (1 - D$lambda)/3 * ASNP[["H10"]]

    objective_val <- lambda_ASN + D$kappa * lambda_ASNP

    if (return_everything) {
      mu_ <- D$mu_vec[["H1"]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[3]), list("TP" = pInf[2], "TC" = pInf[4]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[3]), list("TP" = nInf[2], "TC" = nInf[4]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$A_[["TP1"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP1"]] %*% D$Sigma %*% t(D$A_[["TP1"]]),
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1] + pmvnorm_(
        mean = as.vector(D$A_[["TP12"]] %*% D$mu_vec[["H1"]]),
        sigma = D$A_[["TP12"]] %*% D$Sigma %*% t(D$A_[["TP12"]]),
        lower = c(nInf[[1]][["TP"]], D$b[[2]][["TP"]][["efficacy"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]], pInf[[2]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 2,
        #                       releps = 0)
      )[1]

      D$power <- calc_prob_reject_both("H1", D)
      D$min_conditional_power <- 0


      D$alpha_TP <- D$type_I_error
      D$alpha_TC <- D$type_I_error

      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }

#' Design 1
#'
#' @param x
#' @param D
#' @param return_everything
#'
#' @return
#' @export
#'
#' @examples
objective_single_stage <-
  function(x,
           D,
           return_everything = FALSE) {

    D$nonsequential_futility <- FALSE

    # Add input parameters to Design object
    D$cc <- list()
    D$cc[[1]] <- list()
    D$cc[[1]][["T"]] <- 1
    D$cc[[1]][["P"]] <- x[1]
    D$cc[[1]][["C"]] <- x[2]

    D$b <- list()
    D$b[[2]] <- list()
    D$b[[1]][["TP"]][["efficacy"]] <- qnorm(1 - D$type_I_error)
    D$b[[1]][["TC"]][["efficacy"]] <- qnorm(1 - D$type_I_error)


    # Calculate covariance matrix Sigma
    rho <- list()
    rho[[1]] <- list()
    for (g in c("P", "C")) {
      for (s in 1:(length(D$cc))) {
        rho[[s]][[paste0("T", g)]] <- sqrt(D$sigma[["T"]] / D$cc[[s]][["T"]] + D$sigma[[g]] / D$cc[[s]][[g]])
      }
    }
    D$rho <- rho

    Sigma <- matrix(0, ncol = 2, nrow = 2)
    Sigma[1, ] <-
      c(
        1,
        D$sigma[["T"]] / (D$cc[[1]][["T"]] * rho[[1]][["TP"]] * rho[[1]][["TC"]])
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
          D$mu[["H0"]][["TP"]] / D$rho[[1]][["TP"]],
          D$mu[["H0"]][["TC"]] / D$rho[[1]][["TC"]]
        )
      mu_vec[["H1"]] <-
        sqrt_nT * c(
          D$mu[["H1"]][["TP"]] / D$rho[[1]][["TP"]],
          D$mu[["H1"]][["TC"]] / D$rho[[1]][["TC"]]
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
        algorithm = Miwa(steps = D$maxpts)
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
    D$n[[1]][["P"]] <- nT * D$cc[[1]][["P"]]
    D$n[[1]][["C"]] <- nT * D$cc[[1]][["C"]]

    sqrt_nT <- sqrt(nT)
    mu_vec <- list()
    mu_vec[["H0"]] <-
      sqrt_nT * c(
        D$mu[["H0"]][["TP"]] / D$rho[[1]][["TP"]],
        D$mu[["H0"]][["TC"]] / D$rho[[1]][["TC"]]
      )
    mu_vec[["H1"]] <-
      sqrt_nT * c(
        D$mu[["H1"]][["TP"]] / D$rho[[1]][["TP"]],
        D$mu[["H1"]][["TC"]] / D$rho[[1]][["TC"]]
      )

    D$mu_vec <- mu_vec

    finalStateProbs <- list()
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
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$mvtnorm_alg_maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]

      P[["TP1F"]] <- pmvnorm_(
        mean = as.vector(mu_vec[[hyp]])[1],
        sigma = D$Sigma[1, 1],
        lower = c(nInf[[1]][["TP"]]),
        upper = c(D$b[[1]][["TP"]][["efficacy"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$mvtnorm_alg_maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]

      P[["TP1E_TC1F"]] <- pmvnorm_(
        mean = as.vector(mu_vec[[hyp]]),
        sigma = D$Sigma,
        lower = c(D$b[[1]][["TP"]][["efficacy"]], nInf[[1]][["TC"]]),
        upper = c(pInf[[1]][["TP"]], D$b[[1]][["TC"]][["efficacy"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol / 3,
        #                       releps = 0)
      )[1]
      finalStateProbs[[hyp]] <- P
      ASN[[hyp]] <- sum(unlist(D$n))
      ASNP[[hyp]] <- D$n[[1]][["P"]]
    }

    lambda_ASN <- D$lambda * ASN[["H1"]] + (1 - D$lambda) * ASN[["H0"]]
    lambda_ASNP <- D$lambda * ASNP[["H1"]] + (1 - D$lambda) * ASNP[["H0"]]
    objective_val <- lambda_ASN + D$kappa * lambda_ASNP

    if (return_everything) {
      mu_ <- D$mu_vec[[hyp]]
      pInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = FALSE)
      pInf <- list(list("TP" = pInf[1], "TC" = pInf[2]))
      nInf <- qnorm(D$tol * 1e-2, mean = mu_, lower.tail = TRUE)
      nInf <- list(list("TP" = nInf[1], "TC" = nInf[2]))

      D$x <- x

      D$power_TP <- pmvnorm_(
        mean = as.vector(D$mu_vec[["H1"]][1]),
        sigma = D$Sigma[1, 1],
        lower = c(D$b[[1]][["TP"]][["efficacy"]]),
        upper = c(pInf[[1]][["TP"]]),
        algorithm = Miwa(steps = D$maxpts)
        # algorithm = GenzBretz(maxpts = D$maxpts,
        #                       abseps = D$tol,
        #                       releps = 0)
      )[1]

      D$power <- power
      D$min_conditional_power <- NA_real_

      D$alpha_TP <- D$type_I_error
      D$alpha_TC <- D$type_I_error

      D$finalStateProbs <- finalStateProbs
      D$ASN <- ASN
      D$ASNP <- ASNP
      D$lambda_ASN <- lambda_ASN
      D$lambda_ASNP <- lambda_ASNP
      D$objective_val <- objective_val
      D$maxn <- sum(unlist(D$n))
      return(D)
    } else {
      return(objective_val)
    }
  }
