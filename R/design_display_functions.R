#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.GoldStandardDesign <- function(x, ...){

  sample_sizes_stage1 <- sprintf("Sample sizes (stage 1): T: %i, P: %i, C: %i\n",
                                 x$n[[1]][["T"]], x$n[[1]][["P"]], x$n[[1]][["C"]])
  sample_sizes_stage2 <- sprintf("Sample sizes (stage 2): T: %i, P: %i, C: %i\n",
                                 x$n[[2]][["T"]], x$n[[2]][["P"]], x$n[[2]][["C"]])
  efficacy_boundaries_stage1 <- sprintf("Efficacy boundaries (stage 1): Z_TP_e: %.4f, Z_TC_e: %.4f\n",
                                        x$b[[1]][["TP"]][["efficacy"]], x$b[[1]][["TC"]][["efficacy"]])
  
  futility_boundaries_stage1 <- sprintf("Futility boundaries (stage 1): Z_TP_f: %.4f, Z_TC_f: %.4f\n",
                                        x$b[[1]][["TP"]][["futility"]], x$b[[1]][["TC"]][["futility"]])
  
  efficacy_boundaries_stage2 <- sprintf("Efficacy boundaries (stage 2): Z_TP_e: %.4f, Z_TC_e: %.4f\n",
                                        x$b[[2]][["TP"]][["efficacy"]], x$b[[2]][["TC"]][["efficacy"]])
  max_sample_size <- sprintf("Maximum overall sample size: %i\n", sum(unlist(x$n)))
  ASN_H1 <- sprintf("Expected sample size (H1): %.1f\n", x$ASN$H1)
  ASN_H0 <- sprintf("Expected sample size (H0): %.1f\n", x$ASN$H0)
  cat(sample_sizes_stage1)
  cat(sample_sizes_stage2)
  cat(efficacy_boundaries_stage1)
  cat(futility_boundaries_stage1)
  cat(efficacy_boundaries_stage2)
  cat(max_sample_size)
  cat(ASN_H1)
  cat(ASN_H0)
}

