fr <- function(x, k=2) if (is.character(x)) return(x) else if (is.na(x)) return("") else return(format(round(x, k), nsmall=k))

#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
make_table_2 <- function(D_list){
  table_tib <- tibble(
    `Design`     =     character(),
    `$n_{1, T}$` =     character(),
    `$n_{1, P}$` =     character(),
    `$n_{1, C}$` =     character(),
    `$n_{2, T}$` =     character(),
    `$n_{2, P}$` =     character(),
    `$n_{2, C}$` =     character(),
    `$n_{\\max}$` =    character(),
    `$N_{H_1}$`   =    character(),
    `$CP_{\\min}$`   = character(),
    `$b_{1, TP, f}$` = character(),
    `$b_{1, TP, e}$` = character(),
    `$b_{1, TC, f}$` = character(),
    `$b_{1, TC, e}$` = character(),
    `$b_{2, TP, e}$` = character(),
    `$b_{2, TC, e}$` = character()
    )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        i,
        D$n[[1]][["T"]],
        D$n[[1]][["P"]],
        D$n[[1]][["C"]],
        NA_real_,
        NA_real_,
        NA_real_,
        sum(unlist(D$n)),
        D$ASN[["H1"]],
        NA_real_,
        NA_real_,
        D$b[[1]][["TP"]][["efficacy"]],
        NA_real_,
        D$b[[1]][["TC"]][["efficacy"]],
        NA_real_,
        NA_real_
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(rep(0, 9), rep(2, 7)))
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        i,
        D$cumn[[1]][["T"]],
        D$cumn[[1]][["P"]],
        D$cumn[[1]][["C"]],
        D$cumn[[2]][["T"]],
        D$cumn[[2]][["P"]],
        D$cumn[[2]][["C"]],
        sum(unlist(D$n)),
        D$ASN[["H1"]],
        D$cp_min,
        D$b[[1]][["TP"]][["futility"]],
        D$b[[1]][["TP"]][["efficacy"]],
        if(is.finite(D$b[[1]][["TC"]][["futility"]])) D$b[[1]][["TC"]][["futility"]] else "-\\infty",
        D$b[[1]][["TC"]][["efficacy"]],
        D$b[[2]][["TP"]][["efficacy"]],
        D$b[[2]][["TC"]][["efficacy"]]
        )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(rep(0, 9), rep(2, 7)))
      table_tib <- bind_rows(table_tib, tmp)
    } else{
     stop("Wrong class.")
    }
  }
  return(table_tib)
}

#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
make_table_3 <- function(D_list){
  table_tib <- tibble(
    `$\\lambda$`     =     character(),
    `$n_{1, T}$` =     character(),
    `$n_{1, P}$` =     character(),
    `$n_{1, C}$` =     character(),
    `$n_{2, T}$` =     character(),
    `$n_{2, P}$` =     character(),
    `$n_{2, C}$` =     character(),
    `$n_{\\max}$` =    character(),
    `$N_{H_0}$`   =    character(),
    `$N_{H_1}$`   =    character(),
    `$g_{\\lambda, \\kappa}(D)$` =    character(),
    `$CP_{\\min}$`   = character()
  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        D$lambda,
        D$n[[1]][["T"]],
        D$n[[1]][["P"]],
        D$n[[1]][["C"]],
        NA_real_,
        NA_real_,
        NA_real_,
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$objective_val,
        NA_real_
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, rep(0, 10), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        D$lambda,
        D$cumn[[1]][["T"]],
        D$cumn[[1]][["P"]],
        D$cumn[[1]][["C"]],
        D$cumn[[2]][["T"]],
        D$cumn[[2]][["P"]],
        D$cumn[[2]][["C"]],
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$objective_val,
        D$cp_min
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, rep(0, 10), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else{
      stop("Wrong class.")
    }
  }
  return(table_tib)
}

#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
make_table_4 <- function(D_list){
  table_tib <- tibble(
    `$\\kappa$`     =     character(),
    `$n_{1, T}$` =     character(),
    `$n_{1, P}$` =     character(),
    `$n_{1, C}$` =     character(),
    `$n_{2, T}$` =     character(),
    `$n_{2, P}$` =     character(),
    `$n_{2, C}$` =     character(),
    `$n_{\\max}$` =    character(),
    `$N_{H_0}$`   =    character(),
    `$N_{H_1}$`   =    character(),
    `$N_{H_0}^{P}$`   =    character(),
    `$N_{H_1}^{P}$`   =    character(),
    `$g_{\\lambda, \\kappa}(D)$` =    character(),
    `$CP_{\\min}$`   = character()
  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        "",
        D$n[[1]][["T"]],
        D$n[[1]][["P"]],
        D$n[[1]][["C"]],
        NA_real_,
        NA_real_,
        NA_real_,
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$ASNP[["H0"]],
        D$ASNP[["H1"]],
        D$objective_val,
        NA_real_
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, rep(0, 12), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        D$kappa,
        D$cumn[[1]][["T"]],
        D$cumn[[1]][["P"]],
        D$cumn[[1]][["C"]],
        D$cumn[[2]][["T"]],
        D$cumn[[2]][["P"]],
        D$cumn[[2]][["C"]],
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$ASNP[["H0"]],
        D$ASNP[["H1"]],
        D$objective_val,
        D$cp_min
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, rep(0, 12), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else{
      stop("Wrong class.")
    }
  }
  return(table_tib)
}

#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
make_table_5 <- function(D_list){
  table_tib <- tibble(
    `$\\lambda$`     = character(),
    `$\\kappa$`     =  character(),
    `$n_{1, T}$` =     character(),
    `$n_{1, P}$` =     character(),
    `$n_{1, C}$` =     character(),
    `$n_{2, T}$` =     character(),
    `$n_{2, P}$` =     character(),
    `$n_{2, C}$` =     character(),
    `$n_{\\max}$` =    character(),
    `$N_{H_0}$`   =    character(),
    `$N_{H_1}$`   =    character(),
    `$N_{H_0}^{P}$`   =    character(),
    `$N_{H_1}^{P}$`   =    character(),
    `$g_{\\lambda, \\kappa}(D)$` =    character(),
    `$CP_{\\min}$`   = character()
  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        "",
        D$kappa,
        D$n[[1]][["T"]],
        D$n[[1]][["P"]],
        D$n[[1]][["C"]],
        NA_real_,
        NA_real_,
        NA_real_,
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$ASNP[["H0"]],
        D$ASNP[["H1"]],
        D$objective_val,
        NA_real_
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, 1, rep(0, 12), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        D$lambda,
        D$kappa,
        D$cumn[[1]][["T"]],
        D$cumn[[1]][["P"]],
        D$cumn[[1]][["C"]],
        D$cumn[[2]][["T"]],
        D$cumn[[2]][["P"]],
        D$cumn[[2]][["C"]],
        sum(unlist(D$n)),
        D$ASN[["H0"]],
        D$ASN[["H1"]],
        D$ASNP[["H0"]],
        D$ASNP[["H1"]],
        D$objective_val,
        D$cp_min
      )
      names(tmp) <- names(table_tib)
      tmp <- Map(fr, x=tmp, k=c(1, 1, rep(0, 12), 2))
      table_tib <- bind_rows(table_tib, tmp)
    } else{
      stop("Wrong class.")
    }
  }
  return(table_tib)
}

ccc_wrt_nmax <- function(ccc, maxn, n, singlestage = FALSE){
  newfac <- n[[1]][["T"]] / maxn
  erg <- lapply(ccc, function(x) lapply(x, function(x)x * newfac))
  if (singlestage){
    erg[[2]] <- list("T"=NA_real_, "P"=NA_real_,"C"=NA_real_)
  }
  return(erg)
}

# This should be called TP, typo...
calc_prob_reject_first_stage <- function(hypothesis = "H1", D, groups = "TP"){
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
