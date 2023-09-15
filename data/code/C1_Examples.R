# Example I (Point estimates, confidence intervals, p-values)
set.seed(1)
mean_sirolimus <- 19
sd_sirolimus <- 124
mean_placebo <-  -134
sd_placebo <- 182
gsd <- rpact::getDesignGroupSequential(
  kMax = 2,
  alpha = 0.05, # Yes, the protocol (p. 34) actually states one-sided 5%.
  sided = 1,
  typeOfDesign = "asOF",
  informationRates = c(.4, 1)
)
crit <- qt(pnorm(gsd$criticalValues)[1], df = 38)
nsim <- 100000
erg <- replicate(nsim,
{ # 120 patients with 15-20% attrition rate were assumed in the protocol.
  # Thus, 100 patients were assumed for the simulation.
  x_sirolimus <- rnorm(50, mean = mean_sirolimus, sd = sd_sirolimus)
  x_placebo <- rnorm(50, mean = mean_placebo, sd = sd_placebo)
  if (t.test(x_sirolimus[1:20], x_placebo[1:20])$statistic > crit)
    ret <- mean(x_sirolimus[1:20] - x_placebo[1:20])
  else
    ret <- mean(x_sirolimus - x_placebo)
  ret
})
sprintf(
  "Expected difference: %.4f\nDifference in distributions: %.4f\nBias: %.4f\nRelative bias: %.4f%%",
  mean(erg),
  (mean_sirolimus - mean_placebo),
  mean(erg) - (mean_sirolimus - mean_placebo),
  100*(mean(erg) - (mean_sirolimus - mean_placebo)) / (mean_sirolimus - mean_placebo)
) |> cat()


# Example II (Optimization of gold-standard designs)
n_p <- 243
n_t <- 479
n_c <- 462
mean_p <- -.42
mean_t <- -.61
mean_c <- -.61
sd_p <- .5 # actual observed value was .04
sd_t <- .5 # actual observed value was .03
sd_c <- .5 # actual observed value was .03
noninf_delta <- (mean_p - mean_c) / 2
ci_tp <- mean_t - mean_p + c( -1, 1) * qnorm(.025, lower.tail = FALSE)*sqrt(sd_t^2/n_t + sd_p^2/n_p)
ci_tc <- mean_t - mean_c + c( -1, 1) * qnorm(.025, lower.tail = FALSE)*sqrt(sd_t^2/n_t + sd_c^2/n_c)
sprintf(
"Confidence interval Treatment - Placebo: [%.3f, %.3f]\nConfidence interval Treatment - Control: [%.3f, %.3f]",
ci_tp[1], ci_tp[2], ci_tc[1], ci_tc[2]
) |> cat()

# Create design object according to example
D <- list()
D$n <- list(list(T = n_t, P = n_p, C = n_c))
D$var <- list(T = sd_t^2, P = sd_p^2, C = sd_c^2)
D$mu <- list(H0 = list(TP = 0, TC = 0), H1 = list(TP = mean_p - mean_t, TC = mean_c - mean_t + noninf_delta))
D$stagec <- list(list(T = 1, P = n_p / n_t, C = n_c / n_t))
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
D$mu_wo_nT1 <- list()
for (hyp in c("H0", "H1")){
  D$mu_wo_nT1[[hyp]] <- c(
    D$mu[[hyp]][["TP"]] / D$gamma[[1]][["TP"]],
    D$mu[[hyp]][["TC"]] / D$gamma[[1]][["TC"]]
  )
}
D$b <- list(list(TP = list(efficacy = qnorm(.975)), TC = list(efficacy = qnorm(.975))))
D$mvnorm_algorithm <- mvtnorm::Miwa()
# Calculate power
pwr <- OptimalGoldstandardDesigns:::calc_prob_reject_both_singlestage(D$mu_wo_nT1[["H1"]] * sqrt(D$n[[1]]$T), D)
sprintf("The design from the example has a power of %.0f%%", 100*pwr)

# Optimal single-stage allocation ratios
d1 <- OptimalGoldstandardDesigns::optimize_design_onestage(
  beta = .17, varT = .5^2, varC = .5^2, varP = .5^2,
  Delta = .095, alternative_TP = .19, alternative_TC = 0)
d1
# Optimal Schlömer-Brannath design
d2 <- OptimalGoldstandardDesigns::optimize_design_twostage(
  varT = .5^2, varC = .5^2, varP = .5^2,
  alternative_TP = .19, alternative_TC = 0,
  Delta = .095,
  beta = .17,
  bTP1f = -Inf, bTC1f = -Inf,
  cP1 = d1$stagec[[1]][["P"]], cC1 = d1$stagec[[1]][["C"]],
  cT2 = 1, cP2 = quote(cP1), cC2 = quote(cC1),
  binding_futility = FALSE)
d2
# Design with futility boundaries
d3 <- OptimalGoldstandardDesigns::optimize_design_twostage(
  varT = .5^2, varC = .5^2, varP = .5^2,
  alternative_TP = .19, alternative_TC = 0,
  Delta = .095,
  beta = .17,
  binding_futility = FALSE)
d3
