#' pmvnorm which returns 0 if any lower boundary is larger than
#' any upper boundary
#'
#' @param upper
#' @param lower
#' @param ...
#'
#' @return
#'
#' @examples
#' @importFrom mvtnorm pmvnorm
pmvnorm_ <- function(upper, lower, ...) {
  if (any(lower >= upper)) {
    return(0)
  } else {
    pmvnorm(upper = upper, lower = lower, ...)
  }
}

## Everything below this line is not really relevant for optimal designs
## Some of these functions are used for presentational purposes

emean_to_Z <- function(emean_TP, emean_TC, gamma, nT1) {
  sqrt_nT1 <- sqrt(nT1)
  sqrt_nT1 * c(
    TP = emean_TP / gamma[[1]][["TP"]],
    TC = emean_TC / gamma[[1]][["TC"]]
  )
}

Z_to_emean <- function(Z, gamma, nT1) {
  sqrt_nT1 <- sqrt(nT1)
  list(
    TP = Z[[1]] * gamma[[1]][["TP"]] / sqrt_nT1,
    TC = Z[[2]] * gamma[[1]][["TC"]] / sqrt_nT1
  )
}

conditionalMean <- function(x_a, mu_a, mu_b, Sigma) {
  Sigma_a <- Sigma[seq_len(length(mu_a)), seq_len(length(mu_a))]
  # Sigma_b <- Sigma[(length(mu_a) + 1):(length(mu_a) + length(mu_b)), (length(mu_a) + 1):(length(mu_a) + length(mu_b))]
  Sigma_c <- Sigma[seq_len(length(mu_a)), (length(mu_a) + 1):(length(mu_a) + length(mu_b))]
  return(mu_b + t(Sigma_c) %*% solve(Sigma_a) %*% (x_a - mu_a))
}

conditionalSigma <- function(x_a, mu_a, mu_b, Sigma) {
  Sigma_a <- Sigma[seq_len(length(mu_a)), seq_len(length(mu_a))]
  Sigma_b <- Sigma[(length(mu_a) + 1):(length(mu_a) + length(mu_b)), (length(mu_a) + 1):(length(mu_a) + length(mu_b))]
  Sigma_c <- Sigma[seq_len(length(mu_a)), (length(mu_a) + 1):(length(mu_a) + length(mu_b))]
  return(Sigma_b - t(Sigma_c) %*% solve(Sigma_a) %*% Sigma_c)
}

calc_conditional_power <- function(Z_TP, Z_TC, mu, Sigma, b, nonsequential_futility = FALSE) {
  if (Z_TP <= b[[1]][["TP"]][["futility"]]) {
    return(0)
  } else if (b[[1]][["TP"]][["efficacy"]] < Z_TP) {
    if (Z_TC <= b[[1]][["TC"]][["futility"]]) {
      return(0)
    } else if (b[[1]][["TC"]][["efficacy"]] < Z_TC) {
      return(1)
    } else {
      permutation <- diag(4)[c(1, 3, 2, 4), ]
      cmean <- as.vector(conditionalMean(
        c(Z_TP, Z_TC),
        (permutation %*% mu)[1:2],
        (permutation %*% mu)[3:4],
        permutation %*% Sigma %*% t(permutation)
      ))
      csigma <- conditionalSigma(
        c(Z_TP, Z_TC),
        (permutation %*% mu)[1:2], (permutation %*% mu)[3:4],
        permutation %*% Sigma %*% t(permutation)
      )

      return(pmvnorm_(
        mean = cmean[[2]],
        sigma = csigma[2, 2],
        lower = c(b[[2]][["TC"]][["efficacy"]]),
        upper = c(Inf)
      )[1])
    }
  } else {
    if (isTRUE(nonsequential_futility)) {
      if (Z_TC <= b[[1]][["TC"]][["futility"]]) {
        return(0)
      }
    }
    permutation <- diag(4)[c(1, 3, 2, 4), ]
    cmean <- as.vector(conditionalMean(
      c(Z_TP, Z_TC),
      (permutation %*% mu)[1:2],
      (permutation %*% mu)[3:4], permutation %*% Sigma %*% t(permutation)
    ))
    csigma <- conditionalSigma(
      c(Z_TP, Z_TC), (permutation %*% mu)[1:2],
      (permutation %*% mu)[3:4],
      permutation %*% Sigma %*% t(permutation)
    )
    return(pmvnorm_(
      mean = cmean,
      sigma = csigma,
      lower = c(b[[2]][["TP"]][["efficacy"]], b[[2]][["TC"]][["efficacy"]]),
      upper = c(Inf, Inf)
    )[1])
  }
}

calc_conditional_power_wrt_emean <- function(emean_TP, emean_TC, gamma, nT1, mu, Sigma, b, nonsequential_futility = FALSE) {
  Z <- emean_to_Z(emean_TP, emean_TC, gamma, nT1)
  Z_TP <- Z[["TP"]]
  Z_TC <- Z[["TC"]]
  calc_conditional_power(Z_TP, Z_TC, mu, Sigma, b, nonsequential_futility)
}

calc_conditional_alpha <- function(groups = "TC", Z, mu, Sigma, b, nonsequential_futility = FALSE) {
  stopifnot(groups %in% c("TP", "TC"))

  if (Z <= b[[1]][[groups]][["futility"]]) {
    return(0)
  } else if (b[[1]][[groups]][["efficacy"]] < Z) {
    return(1)
  } else {
    if (groups == "TP") {
      permutation <- diag(4)[c(1, 2, 3, 4), ]
    } else {
      permutation <- diag(4)[c(3, 4, 1, 2), ]
    }
    cmean <- as.vector(conditionalMean(c(Z), (permutation %*% mu)[1], (permutation %*% mu)[2], permutation %*% Sigma %*% t(permutation)))
    csigma <- conditionalSigma(c(Z), (permutation %*% mu)[1], (permutation %*% mu)[2], permutation %*% Sigma %*% t(permutation))

    return(pmvnorm_(
      mean = cmean[[1]],
      sigma = csigma[1, 1],
      lower = c(b[[2]][[groups]][["efficacy"]]),
      upper = c(Inf)
    )[1])
  }
}
