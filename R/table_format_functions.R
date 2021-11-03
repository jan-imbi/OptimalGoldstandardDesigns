ccc_wrt_nmax <- function(ccc, maxn, n, singlestage = FALSE){
  newfac <- n[[1]][["T"]] / maxn
  erg <- lapply(ccc, function(x) lapply(x, function(x)x * newfac))
  if (singlestage){
    erg[[2]] <- list("T"=NA_real_, "P"=NA_real_,"C"=NA_real_)
  }
  return(erg)
}

# This should be called TP, typo...
calc_prob_TC1E <- function(hypothesis = "H1", D, groups = "TP"){
  A_ <- D$A_
  mu_ <- D$mu_vec[[hypothesis]]
  Sigma <- D$Sigma
  b <- D$b
  if (dim(D$Sigma)[1]==2){
    pmvnorm_(
      mean = as.vector(mu_[1]),
      sigma = Sigma[1,1],
      lower = c(b[[1]][[groups]][["efficacy"]]),
      upper = c(Inf),
      algorithm = GenzBretz(maxpts = D$maxpts,
                            abseps = D$tol / 2,
                            releps = 0)
    )[1]
  } else{
    pmvnorm_(
      mean = as.vector(A_[[paste0(groups, "1")]] %*% mu_),
      sigma =  A_[[paste0(groups, "1")]] %*% Sigma %*% t(A_[[paste0(groups, "1")]]),
      lower = c(b[[1]][[groups]][["efficacy"]]),
      upper = c(Inf),
      algorithm = GenzBretz(maxpts = D$maxpts,
                            abseps = D$tol / 2,
                            releps = 0)
    )[1]
  }
}

calc_prob2 <- function(hypothesis = "H1", D, groups = "TP"){
  A_ <- D$A_
  mu_ <- D$mu_vec[[hypothesis]]
  Sigma <- D$Sigma
  b <- D$b
  if (dim(D$Sigma)[1]==2){
    pmvnorm_(
      mean = as.vector(mu_[1]),
      sigma = Sigma[1,1],
      lower = c(b[[1]][[groups]][["efficacy"]]),
      upper = c(Inf),
      algorithm = GenzBretz(maxpts = D$maxpts,
                            abseps = D$tol / 2,
                            releps = 0)
    )[1]
  } else{
    pmvnorm_(
      mean = as.vector(A_[[paste0(groups, "12")]] %*% mu_),
      sigma =  A_[[paste0(groups, "12")]] %*% Sigma %*% t(A_[[paste0(groups, "12")]]),
      lower = c(b[[1]][[groups]][["futility"]], b[[2]][[groups]][["efficacy"]]),
      upper = c(b[[1]][[groups]][["efficacy"]], Inf),
      algorithm = GenzBretz(maxpts = D$maxpts,
                            abseps = D$tol / 2,
                            releps = 0)
    )[1]
  }
}





calc_futility_prob <- function(hypothesis = "H1",
                               groups = "TP",
                               D) {
  if (dim(D$Sigma)[1] == 2) {
    return(0)
  } else{
    A_ <- D$A_
    mu_ <- D$mu_vec[[hypothesis]]
    Sigma <- D$Sigma
    b <- D$b
    return(
      pmvnorm_(
        mean = as.vector(A_[[paste0(groups, "1")]] %*% mu_),
        sigma =  A_[[paste0(groups, "1")]] %*% Sigma %*% t(A_[[paste0(groups, "1")]]),
        lower = c(-Inf),
        upper = c(b[[1]][[groups]][["futility"]]),
        algorithm = GenzBretz(
          maxpts = D$maxpts,
          abseps = D$tol / toladjust,
          releps = 0
        )
      )[1]
    )
  }
}

#
# cumn_wrt_nmax <- function(ccc, maxn, n, singlestage = FALSE){
#   newfac <- n[[1]][["T"]] / maxn
#   erg <- lapply(ccc, function(x) lapply(x, function(x)x * newfac))
#   if (ss){
#     erg[[2]] <- list("T"=NA_real_, "P"=NA_real_,"C"=NA_real_)
#   }
#   return(erg)
# }
