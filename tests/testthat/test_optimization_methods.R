a <- optimize_design_twostage(nloptr_opts = list(
  algorithm = "NLOPT_LN_SBPLX",
  xtol_rel = .1, # inner_tol_objective / 9
  maxeval = 20
))


b <- optimize_design_onestage()
