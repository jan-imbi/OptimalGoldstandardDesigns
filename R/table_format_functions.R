fr <- function(x, k=2){
  format_single <- function(x){
    if (is.character(x))
      return(x)
    else if (is.na(x))
      return("")
    else if (is.infinite(x))
      if (x > 0){
        return("$\\infty$")
      } else{
        return("$-\\infty$")
      }
    else return(format(round(x, k), nsmall=k))
  }
  sapply(x, format_single)
}
#' @importFrom scales percent
ri <- function(d1, d2){
  scales::percent(
    1 - d2$ASN$H1 /
      d1$ASN$H1, accuracy=.1)
}
#' @importFrom stats pnorm
pTP1E <- function(D){
  mu_ <- D$mu_vec[["H1"]]
  b <- D$b
  scales::percent(
  pnorm(b[[1]][["TP"]][["efficacy"]],
        mean = mu_[[1]], sd=1, lower.tail = FALSE),
  accuracy=.1)
}
#' @importFrom stats pnorm
pTP1E_TC1E <- function(D){
  scales::percent(
    D$final_state_probs$H1$TP1E_TC1E,
    accuracy=.1)
}

#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
make_table_2 <- function(D_list){
  # Map(fr, x=tmp, k=c(rep(0, 9), rep(2, 7)))
  table_tib <- tibble(
    `Design`         = numeric(),
    `$n_{1, T}$`     = numeric(),
    `$n_{1, P}$`     = numeric(),
    `$n_{1, C}$`     = numeric(),
    `$n_{2, T}$`     = numeric(),
    `$n_{2, P}$`     = numeric(),
    `$n_{2, C}$`     = numeric(),
    `$n_{\\max}$`    = numeric(),
    `$N_{H_1}$`      = numeric(),
    `$CP_{\\min}$`   = numeric(),
    `$b_{1, TP, f}$` = numeric(),
    `$b_{1, TP, e}$` = numeric(),
    `$b_{1, TC, f}$` = numeric(),
    `$b_{1, TC, e}$` = numeric(),
    `$b_{2, TP, e}$` = numeric(),
    `$b_{2, TC, e}$` = numeric()
    )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
         `Design`         = if(!is.null(D$nr)) D$nr else NA_real_,
         `$n_{1, T}$`     = D$n[[1]][["T"]],
         `$n_{1, P}$`     = D$n[[1]][["P"]],
         `$n_{1, C}$`     = D$n[[1]][["C"]],
         `$n_{2, T}$`     = NA_real_,
         `$n_{2, P}$`     = NA_real_,
         `$n_{2, C}$`     = NA_real_,
         `$n_{\\max}$`    = sum(unlist(D$n)),
         `$N_{H_1}$`      = D$ASN[["H1"]],
         `$CP_{\\min}$`   = NA_real_,
         `$b_{1, TP, f}$` = NA_real_,
         `$b_{1, TP, e}$` = D$b[[1]][["TP"]][["efficacy"]],
         `$b_{1, TC, f}$` = NA_real_,
         `$b_{1, TC, e}$` = D$b[[1]][["TC"]][["efficacy"]],
         `$b_{2, TP, e}$` = NA_real_,
         `$b_{2, TC, e}$` = NA_real_
      )
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
         `Design`         = if(!is.null(D$nr)) D$nr else NA_real_,
         `$n_{1, T}$`     = D$cumn[[1]][["T"]],
         `$n_{1, P}$`     = D$cumn[[1]][["P"]],
         `$n_{1, C}$`     = D$cumn[[1]][["C"]],
         `$n_{2, T}$`     = D$cumn[[2]][["T"]],
         `$n_{2, P}$`     = D$cumn[[2]][["P"]],
         `$n_{2, C}$`     = D$cumn[[2]][["C"]],
         `$n_{\\max}$`    = sum(unlist(D$n)),
         `$N_{H_1}$`      = D$ASN[["H1"]],
         `$CP_{\\min}$`   = D$cp_min,
         `$b_{1, TP, f}$` = D$b[[1]][["TP"]][["futility"]],
         `$b_{1, TP, e}$` = D$b[[1]][["TP"]][["efficacy"]],
         `$b_{1, TC, f}$` = D$b[[1]][["TC"]][["futility"]],
         `$b_{1, TC, e}$` = D$b[[1]][["TC"]][["efficacy"]],
         `$b_{2, TP, e}$` = D$b[[2]][["TP"]][["efficacy"]],
         `$b_{2, TC, e}$` = D$b[[2]][["TC"]][["efficacy"]]
        )
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
  # Map(fr, x=tmp, k=c(1, rep(0, 10), 2))
  table_tib <- tibble(
    `$\\lambda$`                 = numeric(),
    `$n_{1, T}$`                 = numeric(),
    `$n_{1, P}$`                 = numeric(),
    `$n_{1, C}$`                 = numeric(),
    `$n_{2, T}$`                 = numeric(),
    `$n_{2, P}$`                 = numeric(),
    `$n_{2, C}$`                 = numeric(),
    `$n_{\\max}$`                = numeric(),
    `$N_{H_0}$`                  = numeric(),
    `$N_{H_1}$`                  = numeric(),
    `$g_{\\lambda, \\kappa}(D)$` = numeric(),
    `$CP_{\\min}$`               = numeric(),
    `$b_{1, TP, f}$`             = numeric(),
    `$b_{1, TC, f}$`             = numeric()

  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
       `$\\lambda$`                 =  D$lambda,
       `$n_{1, T}$`                 =  D$n[[1]][["T"]],
       `$n_{1, P}$`                 =  D$n[[1]][["P"]],
       `$n_{1, C}$`                 =  D$n[[1]][["C"]],
       `$n_{2, T}$`                 =  NA_real_,
       `$n_{2, P}$`                 =  NA_real_,
       `$n_{2, C}$`                 =  NA_real_,
       `$n_{\\max}$`                =  sum(unlist(D$n)),
       `$N_{H_0}$`                  =  D$ASN[["H0"]],
       `$N_{H_1}$`                  =  D$ASN[["H1"]],
       `$g_{\\lambda, \\kappa}(D)$` =  D$objective_val,
       `$CP_{\\min}$`               =  NA_real_,
       `$b_{1, TP, f}$`             =  D$b[[1]][["TP"]][["futility"]],
       `$b_{1, TC, f}$`             =  D$b[[1]][["TC"]][["futility"]]
      )
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        `$\\lambda$`                 =  D$lambda,
        `$n_{1, T}$`                 =  D$cumn[[1]][["T"]],
        `$n_{1, P}$`                 =  D$cumn[[1]][["P"]],
        `$n_{1, C}$`                 =  D$cumn[[1]][["C"]],
        `$n_{2, T}$`                 =  D$cumn[[2]][["T"]],
        `$n_{2, P}$`                 =  D$cumn[[2]][["P"]],
        `$n_{2, C}$`                 =  D$cumn[[2]][["C"]],
        `$n_{\\max}$`                =  sum(unlist(D$n)),
        `$N_{H_0}$`                  =  D$ASN[["H0"]],
        `$N_{H_1}$`                  =  D$ASN[["H1"]],
        `$g_{\\lambda, \\kappa}(D)$` =  D$objective_val,
        `$CP_{\\min}$`               =  D$cp_min,
        `$b_{1, TP, f}$`             =  D$b[[1]][["TP"]][["futility"]],
        `$b_{1, TC, f}$`             =  D$b[[1]][["TC"]][["futility"]]
      )
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
  # Map(fr, x=tmp, k=c(1, rep(0, 12), 2))
  table_tib <- tibble(
    `$\\kappa$`                  = numeric(),
    `$n_{1, T}$`                 = numeric(),
    `$n_{1, P}$`                 = numeric(),
    `$n_{1, C}$`                 = numeric(),
    `$n_{2, T}$`                 = numeric(),
    `$n_{2, P}$`                 = numeric(),
    `$n_{2, C}$`                 = numeric(),
    `$n_{\\max}$`                = numeric(),
    `$N_{H_0}$`                  = numeric(),
    `$N_{H_1}$`                  = numeric(),
    `$N_{H_0}^{P}$`              = numeric(),
    `$N_{H_1}^{P}$`              = numeric(),
    `$g_{\\lambda, \\kappa}(D)$` = numeric(),
    `$CP_{\\min}$`               = numeric()
  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        `$\\kappa$`                  =  "",
        `$n_{1, T}$`                 =  D$n[[1]][["T"]],
        `$n_{1, P}$`                 =  D$n[[1]][["P"]],
        `$n_{1, C}$`                 =  D$n[[1]][["C"]],
        `$n_{2, T}$`                 =  NA_real_,
        `$n_{2, P}$`                 =  NA_real_,
        `$n_{2, C}$`                 =  NA_real_,
        `$n_{\\max}$`                =  sum(unlist(D$n)),
        `$N_{H_0}$`                  =  D$ASN[["H0"]],
        `$N_{H_1}$`                  =  D$ASN[["H1"]],
        `$N_{H_0}^{P}$`              =  D$ASNP[["H0"]],
        `$N_{H_1}^{P}$`              =  D$ASNP[["H1"]],
        `$g_{\\lambda, \\kappa}(D)$` =  D$objective_val,
        `$CP_{\\min}$`               =  NA_real_
      )
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        `$\\kappa$`                  =  D$kappa,
        `$n_{1, T}$`                 =  D$cumn[[1]][["T"]],
        `$n_{1, P}$`                 =  D$cumn[[1]][["P"]],
        `$n_{1, C}$`                 =  D$cumn[[1]][["C"]],
        `$n_{2, T}$`                 =  D$cumn[[2]][["T"]],
        `$n_{2, P}$`                 =  D$cumn[[2]][["P"]],
        `$n_{2, C}$`                 =  D$cumn[[2]][["C"]],
        `$n_{\\max}$`                =  sum(unlist(D$n)),
        `$N_{H_0}$`                  =  D$ASN[["H0"]],
        `$N_{H_1}$`                  =  D$ASN[["H1"]],
        `$N_{H_0}^{P}$`              =  D$ASNP[["H0"]],
        `$N_{H_1}^{P}$`              =  D$ASNP[["H1"]],
        `$g_{\\lambda, \\kappa}(D)$` =  D$objective_val,
        `$CP_{\\min}$`               =  D$cp_min
      )
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
  # Map(fr, x=tmp, k=c(1, 1, rep(0, 12), 2))
  table_tib <- tibble(
    `$\\lambda$`                 = numeric(),
    `$\\kappa$`                  = numeric(),
    `$n_{1, T}$`                 = numeric(),
    `$n_{1, P}$`                 = numeric(),
    `$n_{1, C}$`                 = numeric(),
    `$n_{2, T}$`                 = numeric(),
    `$n_{2, P}$`                 = numeric(),
    `$n_{2, C}$`                 = numeric(),
    `$n_{\\max}$              `  = numeric(),
    `$N_{H_0}$`                  = numeric(),
    `$N_{H_1}$`                  = numeric(),
    `$N_{H_0}^{P}$`              = numeric(),
    `$N_{H_1}^{P}$`              = numeric(),
    `$g_{\\lambda, \\kappa}(D)$` = numeric(),
    `$CP_{\\min}$`               = numeric()
  )
  for (i in seq_along(D_list)){
    D <- D_list[[i]]
    if ("OneStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
       `$\\lambda$`                 =   "",
       `$\\kappa$`                  =   D$kappa,
       `$n_{1, T}$`                 =   D$n[[1]][["T"]],
       `$n_{1, P}$`                 =   D$n[[1]][["P"]],
       `$n_{1, C}$`                 =   D$n[[1]][["C"]],
       `$n_{2, T}$`                 =   NA_real_,
       `$n_{2, P}$`                 =   NA_real_,
       `$n_{2, C}$`                 =   NA_real_,
       `$n_{\\max}$              `  =   sum(unlist(D$n)),
       `$N_{H_0}$`                  =   D$ASN[["H0"]],
       `$N_{H_1}$`                  =   D$ASN[["H1"]],
       `$N_{H_0}^{P}$`              =   D$ASNP[["H0"]],
       `$N_{H_1}^{P}$`              =   D$ASNP[["H1"]],
       `$g_{\\lambda, \\kappa}(D)$` =   D$objective_val,
       `$CP_{\\min}$`               =   NA_real_
      )
      table_tib <- bind_rows(table_tib, tmp)
    } else if ("TwoStageGoldStandardDesign" %in% class(D)){
      tmp <- tibble(
        `$\\lambda$`                 =  D$lambda,
        `$\\kappa$`                  =  D$kappa,
        `$n_{1, T}$`                 =  D$cumn[[1]][["T"]],
        `$n_{1, P}$`                 =  D$cumn[[1]][["P"]],
        `$n_{1, C}$`                 =  D$cumn[[1]][["C"]],
        `$n_{2, T}$`                 =  D$cumn[[2]][["T"]],
        `$n_{2, P}$`                 =  D$cumn[[2]][["P"]],
        `$n_{2, C}$`                 =  D$cumn[[2]][["C"]],
        `$n_{\\max}$              `  =  sum(unlist(D$n)),
        `$N_{H_0}$`                  =  D$ASN[["H0"]],
        `$N_{H_1}$`                  =  D$ASN[["H1"]],
        `$N_{H_0}^{P}$`              =  D$ASNP[["H0"]],
        `$N_{H_1}^{P}$`              =  D$ASNP[["H1"]],
        `$g_{\\lambda, \\kappa}(D)$` =  D$objective_val,
        `$CP_{\\min}$`               =  D$cp_min
      )
      table_tib <- bind_rows(table_tib, tmp)
    } else{
      stop("Wrong class.")
    }
  }
  return(table_tib)
}
