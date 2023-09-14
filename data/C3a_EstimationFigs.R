# Author: Jan Meis
# Source: https://github.com/jan-imbi/adestr/

# If you did not run the calculations.R script beforehand, you will need to get
# the data from https://github.com/jan-imbi/adestr/tree/master/data
# and save it in the "./data/" directory

# Load all data from the /data/ directory
files <- list.files("data") |>
  grepl(".*\\.(rds|RDS)", x=_) |>
  subset(list.files("data"), subset = _) |>
  (\(x) paste0("data/", x))()
for (f in files){
  assign(stringr::str_extract(f, "(?<=data/).*(?=\\.(rds|RDS))"), readRDS(f))
}
dir.create("figures")

# All packages and helper functions required to produce the plots
library(adestr)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)
library(Cairo)
sanitize_names <- function(x) {
  if (is.data.frame(x)) {
    xx <- x
    x <- names(xx)
  }
  x <- stringr::str_to_sentence(x)
  x[x=="Mse"] <- "MSE"
  x[x=="Mu"] <- "mu"
  x[x=="Sigma"] <- "sigma"
  if (is.data.frame(xx)) {
    names(xx) <- x
    return(xx)
  } else{
    return(x)
  }
}
sanitize_vars <- function(x) {
  xest <- x$Estimator
  xest[xest=="Bias reduced MLE (iterations=1)"] <- "Bias-reduced MLE (1 iteration)"
  xest[xest=="Median unbiased (LR test ordering)"] <- "Median-unbiased (LR test ordering)"
  xest[xest=="Median unbiased (MLE ordering)"] <- "Median-unbiased (MLE ordering)"
  xest[xest=="Median unbiased (Score test ordering)"] <- "Median-unbiased (Score test ordering)"
  xest[xest=="Median unbiased (SWCF ordering)"] <- "Median-unbiased (SWCF ordering)"
  xest[xest=="LR test ordering"] <- "LR test ordering CI"
  xest[xest=="MLE ordering"] <- "MLE ordering CI"
  xest[xest=="Score test ordering"] <- "Score test ordering CI"
  xest[xest=="SWCF ordering"] <- "SWCF ordering CI"
  xest[xest==""]
  x$Estimator <- xest
  x
}

### Figure 4 ###
# Design characteristics of adaptive vs. group-sequential design (normal distribution)
f_designs <- adestr:::plot_design(list(designad, designgs), data_distribution = Normal(two_armed = FALSE))
f_designs

ggsave(filename = "figures/f_designs.pdf",
       plot = f_designs,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))


### Figure 5 ###
# Design characteristics of adaptive vs. group-sequential design (t-distribution)
f_designs_t <- adestr:::plot_design(list(designad_t, designgs_t), data_distribution = Student(two_armed = TRUE))
f_designs_t

ggsave(filename = "figures/f_designs_t.pdf",
       plot = f_designs_t,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))


### Figure 7 ###
# Probability of continuation
df_continuation <- tibble::tibble(
  mu = c(muvec, muvec),
  `Type of design` = rep(c("Adaptive", "Group-sequential"), each = length(muvec)),
  `Probability of continuation` = c(
    pnorm(designad@c1e, mean=muvec * sqrt(designad@n1)) - pnorm(designad@c1f, mean=muvec* sqrt(designad@n1)),
    pnorm(designgs@c1e, mean=muvec * sqrt(designgs@n1)) - pnorm(designgs@c1f, mean=muvec* sqrt(designgs@n1))
  )
)
f_est_cont <-ggplot(df_continuation, mapping = aes(x = mu, y = `Probability of continuation`, col = `Type of design`)) +
  geom_line(linewidth=1) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, .6, .1), limits = c(0, .6)) +
  scale_x_continuous(TeX("$\\mu$")) +
  geom_vline(aes(xintercept = designad@c1f/sqrt(designad@n1))) +
  geom_vline(aes(xintercept = designad@c1e/sqrt(designad@n1))) +
  theme_pubr()
f_est_cont

ggsave(filename = "figures/f_est_cont.pdf",
       plot = f_est_cont,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))


### Figure 8 ###
# Distribution of the sample mean (normal, one-armed case)
f_smean <- adestr:::plot_sample_mean(Normal(two_armed = FALSE), designad, mu = seq(-.3, .8, .1), sigma = 1,
                                     combine_components = TRUE, exact = FALSE) + theme(text = element_text(size=9))
f_smean
saveRDS(f_smean, "data/f_smean.rds")
ggsave(filename = "figures/f_smean.pdf",
       plot = f_smean,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))

### Figures 9, 10, 11, 12 ###
# Estimator characteristics
tab_est_filter <- est_results |>
  dplyr::filter(!(estimator%in% c("MedianUnbiasedNeymanPearsonOrdering(mu0=0, mu1=0.4)"))) |>
  sanitize_names() |>
  sanitize_vars()
f_est_bias <- ggplot(data = tab_est_filter, aes(x = mu, y=Bias, col = Estimator)) +
  geom_line(linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_var <- ggplot() +
  geom_line(data = tab_est_filter, aes(x = mu, y=Variance, col = Estimator), linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_mse <- ggplot() +
  geom_line(data = tab_est_filter, aes(x = mu, y=MSE, col = Estimator), linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_overest <- ggplot(data = tab_est_filter, aes(x = mu, y=`Probability of overestimation`, col = Estimator)) +
  geom_line(linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3)) +
  scale_y_continuous(limits = c(.35, .65), breaks = seq(.35, .65, .05))
f_est_bias
f_est_var
f_est_mse
f_est_overest

ggsave(filename = "figures/f_est_bias.pdf",
       plot = f_est_bias,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))
ggsave(filename = "figures/f_est_var.pdf",
       plot = f_est_var,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))
ggsave(filename = "figures/f_est_mse.pdf",
       plot = f_est_mse,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))
ggsave(filename = "figures/f_est_overest.pdf",
       plot = f_est_overest,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))

### Figures 13, 14, 15 ###
# Confidence interval characteristics
tab_ci_filter <- ci_results |>
  dplyr::filter(!(estimator%in% c(
    "Neyman-Pearson test ordering (mu0=0, mu1=0.4)"
  ))) |>
  sanitize_names() |>
  sanitize_vars()
f_ci_cov <- ggplot(data = tab_ci_filter, aes(x = mu, y=Coverage, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_width <- ggplot(data = tab_ci_filter, aes(x = mu, y=Width, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  scale_y_continuous("Expected width") +
  theme_pubr()
f_ci_test <- ggplot(data = tab_ci_filter, aes(x = mu, y=`Agreement with test decision`, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_cov
f_ci_width
f_ci_test

ggsave(filename = "figures/f_ci_cov.pdf",
       plot = f_ci_cov,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))
ggsave(filename = "figures/f_ci_width.pdf",
       plot = f_ci_width,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))
ggsave(filename = "figures/f_ci_test.pdf",
       plot = f_ci_test,
       device = "pdf",
       width = unit(7, "cm"), height = unit(4, "cm"))


### Figures 16, 17 ###
# P-values / stage-wise ordering representation for an adaptive design
f_p_sw <- adestr:::plot_p(StagewiseCombinationFunctionOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[1], subdivisions = 200L)
f_p_ml <- adestr:::plot_p(MLEOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[2], subdivisions = 200L)
f_p_lr <- adestr:::plot_p(LikelihoodRatioOrderingPValue(), Normal(two_armed = FALSE),  designad, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)
f_p_st <- adestr:::plot_p(ScoreTestOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[4], subdivisions = 200L)
f_p_np <- adestr:::plot_p(NeymanPearsonOrderingPValue(), Normal(two_armed = FALSE), designad, -1, 1, boundary_color = scales::hue_pal()(5)[5], subdivisions = 200L)
f_p_reject <- adestr:::plot_rejection_regions(
  list(StagewiseCombinationFunctionOrderingPValue(), MLEOrderingPValue(),  LikelihoodRatioOrderingPValue(), ScoreTestOrderingPValue(), NeymanPearsonOrderingPValue()),
  Normal(two_armed = FALSE), designad, 0, 1, subdivisions = 200L)
f_p <- ggarrange(f_p_sw, f_p_ml, f_p_lr, f_p_st, f_p_np, f_p_reject,  ncol = 2, nrow = 3,
                 legend.grob = gridExtra::gtable_rbind(get_legend(f_p_reject), get_legend(f_p_sw)))
Cairo::CairoPNG("figures/f_p.png", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
dev.off()

# P-values / stage-wise ordering representation for a group-sequential design
f_p_sw <- adestr:::plot_p(StagewiseCombinationFunctionOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[1], subdivisions = 200L)
f_p_ml <- adestr:::plot_p(MLEOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[2], subdivisions = 200L)
f_p_lr <- adestr:::plot_p(LikelihoodRatioOrderingPValue(), Normal(two_armed = FALSE),  designgs, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)
f_p_st <- adestr:::plot_p(ScoreTestOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[4], subdivisions = 200L)
f_p_np <- adestr:::plot_p(NeymanPearsonOrderingPValue(), Normal(two_armed = FALSE), designgs, -1, 1, boundary_color = scales::hue_pal()(5)[5], subdivisions = 200L)
f_p_reject <- adestr:::plot_rejection_regions(
  list(StagewiseCombinationFunctionOrderingPValue(), MLEOrderingPValue(),  LikelihoodRatioOrderingPValue(), ScoreTestOrderingPValue(), NeymanPearsonOrderingPValue()),
  Normal(two_armed = FALSE), designgs, 0, 1, subdivisions = 200L)
f_p <- ggarrange(f_p_sw, f_p_ml, f_p_lr, f_p_st, f_p_np, f_p_reject,  ncol = 2, nrow = 3,
                 legend.grob = gridExtra::gtable_rbind(get_legend(f_p_reject), get_legend(f_p_sw)))
Cairo::CairoPNG("figures/f_p_gs.png", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
dev.off()


### Figures 18, 19 ###
# Estimator characteristics (t-distribution two-armed case)
tab_est_filter_t <- est_results_t |>
  dplyr::filter(!(estimator%in% c("MedianUnbiasedNeymanPearsonOrdering(mu0=0, mu1=0.4)"))) |>
  sanitize_names() |>
  sanitize_vars()
f_est_mse_t <- ggplot() +
  geom_line(data = tab_est_filter_t, aes(x = mu, y=MSE, col = Estimator), linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_bias_t <- ggplot(data = tab_est_filter_t, aes(x = mu, y=Bias, col = Estimator)) +
  geom_line(linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_t <- ggarrange(f_est_bias_t, f_est_mse_t, common.legend = TRUE)
f_est_t
ggsave(filename = "figures/f_est_t.pdf",
       plot = f_est_t,
       device = "pdf",
       width = unit(7, "cm"), height = unit(3.5, "cm"))

# Confidence interval characteristics
tab_ci_filter_t <- ci_results_t |>
  dplyr::filter(!(estimator%in% c(
    "Neyman-Pearson test ordering (mu0=0, mu1=0.4)"
  ))) |>
  sanitize_names() |>
  dplyr::filter(-0.74 < mu & mu < 1.33)
f_ci_cov_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=Coverage, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_width_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=Width, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  scale_y_continuous("Expected width") +
  theme_pubr()
f_ci_test_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=`Agreement with test decision`, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_t <- ggarrange(f_ci_cov_t, f_ci_width_t, f_ci_test_t, common.legend = TRUE, nrow=3)
f_ci_t
ggsave(filename = "figures/f_ci_t.pdf",
       plot = f_ci_t,
       device = "pdf",
       width = unit(7, "cm"), height = unit(9, "cm"))
