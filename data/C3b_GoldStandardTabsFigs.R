# Instructions (recommended):
# git clone https://github.com/jan-imbi/OptimalGoldstandardDesigns
# open the project and run the script from there.
#
# Instructions (alternative):
# construct the following folder structure:
# root
# |-- R
# |-- data
# |-- |-- code
#
# Put reproduce_everything.R in root/data/code, and all other R scripts in root/R.
# Open a project in the root directory and run reproduce_everything.R

# Load all data from the /data/ directory
files <- list.files("data") |>
  grepl(".*\\.rds", x=_) |>
  subset(list.files("data"), subset = _) |>
  (\(x) paste0("data/", x))()
for (f in files){
  assign(stringr::str_extract(f, "(?<=data/).*(?=\\.rds)"), readRDS(f))
}
dir.create("figures")

library(OptimalGoldstandardDesigns)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(latex2exp)
library(ggpubr)

### Separate the data into tables for the 5 Subchapters ###
tab2 <- goldstandard_results |>
  filter(Design %in% c(1,2,3) | (Design==4) & `$\\kappa$` == 0 | (Design==5) & `$\\kappa$` == 0 & `$\\lambda$` == 1 & `$\\eta$` == 0) |>
  arrange(desc(`$\\beta$`), Design)
tab3 <- goldstandard_results |>
  filter((Design==5) & `$\\kappa$` == 0 & `$\\eta$` == 0 & `$\\lambda$` > 0.15) |>
  arrange(desc(`$\\beta$`), desc(`$\\lambda$`))
tab4 <- goldstandard_results |>
  filter((Design==5) & `$\\lambda$` == 1 & `$\\eta$` == 0) |>
  arrange(desc(`$\\beta$`), `$\\kappa$`)
tab5 <- goldstandard_results |>
  filter((Design==5) & `$\\lambda$` == 1 & `$\\kappa$` == 0) |>
  arrange(desc(`$\\beta$`), `$\\eta$`)
tab6 <- goldstandard_results |>
  filter(Design == 2 |  `$\\lambda$` == .9 & `$\\eta$` == 0.2 | `$\\lambda$` == .95 & `$\\eta$` == 0.05 | `$\\lambda$` == .96 & `$\\eta$` == 0.06) |>
  arrange(desc(`$\\beta$`), `$\\eta$`)


### Numbers in the text in Chapter 3.2.1 (Comparison of the five designs) ###
(tab2$`$N_{H_1}$`[1:5]/ tab2$`$N_{H_1}$`[2]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()
(tab2$`$N_{H_1}$`[6:10]/ tab2$`$N_{H_1}$`[7]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()
(tab2$`$N_{H_0}$`[1:5]/ tab2$`$N_{H_0}$`[2]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()
(tab2$`$N_{H_0}$`[6:10]/ tab2$`$N_{H_0}$`[7]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()

### Relativ sample sizes in chapter "Combining hypothesis weighting with placebo group
### and maximum sample size penalties"
(tab6$`$N_{H_0}$`[2:3]/ tab6$`$N_{H_0}$`[1]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()
(tab6$`$N_{H_0}$`[5:6]/ tab6$`$N_{H_0}$`[4]*100) |>
  OptimalGoldstandardDesigns:::fr(0) |>
  paste0("%, ") |>
  cat()

### Tables 3.2.1
tab2_n <- tab2 |>
  select(Design, `$\\beta$`,  `$n_{1, T}$`, `$n_{1, P}$`, `$n_{1, C}$`, `$n_{2, T}$`, `$n_{2, P}$`, `$n_{2, C}$`)
tab2_b <- tab2 |>
  select(Design, `$\\beta$`,
         `$b_{1, TP, f}$`, `$b_{1, TP, e}$`, `$b_{1, TC, f}$`, `$b_{1, TC, e}$`, `$b_{2, TP, e}$`, `$b_{2, TC, e}$`)
tab2_res <- tab2 |>
  select(Design, `$\\beta$`,
         `$n_{\\max}$`, `$N_{H_0}$`, `$N_{H_1}$`, `$N_{H_0}^{P}$`, `$N_{H_1}^{P}$`,
         `$CP_{\\min}$`)

tn <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab2_n, k = c(0, 1, rep(0, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_2_n} Cumulative sample sizes of the different designs."
)
tb <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab2_b, k = c(0, 1, rep(2, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "",  "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_2_b} Z-scale boundaries of the different designs."
)
tr <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab2_res, k = c(0, 1, 0, 0, 0, 0, 0, 2))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_2_res} Operating characteristics of the different designs."
)
cat(tn, tb, tr, sep = "\n\n")


### Table 3.2.2
tab3_n <- tab3 |>
  select(`$\\beta$`, `$\\lambda$`, `$n_{1, T}$`, `$n_{1, P}$`, `$n_{1, C}$`, `$n_{2, T}$`, `$n_{2, P}$`, `$n_{2, C}$`)
tab3_b <- tab3 |>
  select(`$\\beta$`, `$\\lambda$`,
         `$b_{1, TP, f}$`, `$b_{1, TP, e}$`, `$b_{1, TC, f}$`, `$b_{1, TC, e}$`, `$b_{2, TP, e}$`, `$b_{2, TC, e}$`)
tab3_res <- tab3 |>
  select(`$\\beta$`, `$\\lambda$`,
         `$n_{\\max}$`, `$N_{H_0}$`, `$N_{H_1}$`, `$N_{H_0}^{P}$`, `$N_{H_1}^{P}$`,
         `$CP_{\\min}$`)

tn <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab3_n, k = c(1, 1, rep(0, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "",  "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_3_n} Cumulative sample sizes for different values of $\\lambda$."
)
tb <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab3_b, k = c(1, 1, rep(2, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "", "",  "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_3_b} Z-scale boundaries for different values of $\\lambda$."
)
tr <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab3_res, k = c(1, 1, 0, 0, 0, 0, 0, 2))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "",  "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_3_res} Operating characteristics for different values of $\\lambda$."
)
cat(tn, tb, tr, sep = "\n\n")



### Table 3.2.3
tab4_n <- tab4 |>
  select(`$\\beta$`, `$\\kappa$`, `$n_{1, T}$`, `$n_{1, P}$`, `$n_{1, C}$`, `$n_{2, T}$`, `$n_{2, P}$`, `$n_{2, C}$`)
tab4_b <- tab4 |>
  select(`$\\beta$`, `$\\kappa$`,
         `$b_{1, TP, f}$`, `$b_{1, TP, e}$`, `$b_{1, TC, f}$`, `$b_{1, TC, e}$`, `$b_{2, TP, e}$`, `$b_{2, TC, e}$`)
tab4_res <- tab4 |>
  select(`$\\beta$`, `$\\kappa$`,
         `$n_{\\max}$`, `$N_{H_0}$`, `$N_{H_1}$`, `$N_{H_0}^{P}$`, `$N_{H_1}^{P}$`,
         `$CP_{\\min}$`)
tn <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab4_n, k = c(1, 1, rep(0, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_4_n} Cumulative sample sizes for different values of $\\kappa$."
)
tb <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab4_b, k = c(1, 1, rep(2, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_4_b} Z-scale boundaries for different values of $\\kappa$."
)
tr <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab4_res, k = c(1, 1, 0, 0, 0, 0, 0, 2))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_4_res} Operating characteristics for different values of $\\kappa$."
)
cat(tn, tb, tr, sep = "\n\n")

## Tables 3.2.4
tab5_n <- tab5 |>
  select(`$\\beta$`, `$\\eta$`, `$n_{1, T}$`, `$n_{1, P}$`, `$n_{1, C}$`, `$n_{2, T}$`, `$n_{2, P}$`, `$n_{2, C}$`)
tab5_b <- tab5 |>
  select(`$\\beta$`, `$\\eta$`,
         `$b_{1, TP, f}$`, `$b_{1, TP, e}$`, `$b_{1, TC, f}$`, `$b_{1, TC, e}$`, `$b_{2, TP, e}$`, `$b_{2, TC, e}$`)
tab5_res <- tab5 |>
  select(`$\\beta$`, `$\\eta$`,
         `$n_{\\max}$`, `$N_{H_0}$`, `$N_{H_1}$`, `$N_{H_0}^{P}$`, `$N_{H_1}^{P}$`,
         `$CP_{\\min}$`)
tn <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab5_n, k = c(1, 1, rep(0, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_5_n} Cumulative sample sizes for different values of $\\eta$."
)
tb <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab5_b, k = c(1, 1, rep(2, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_5_b} Z-scale boundaries for different values of $\\eta$."
)
tr <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab5_res, k = c(1,  1, 0, 0, 0, 0, 0, 2))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "", "", "", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_5_res} Operating characteristics for different values of $\\eta$."
)
cat(tn, tb, tr, sep = "\n\n")


### Tables 3.2.5
tab6_n <- tab6 |>
  select(Design, `$\\beta$`, `$n_{1, T}$`, `$n_{1, P}$`, `$n_{1, C}$`, `$n_{2, T}$`, `$n_{2, P}$`, `$n_{2, C}$`)
tab6_b <- tab6 |>
  select(Design, `$\\beta$`,
         `$b_{1, TP, f}$`, `$b_{1, TP, e}$`, `$b_{1, TC, f}$`, `$b_{1, TC, e}$`, `$b_{2, TP, e}$`, `$b_{2, TC, e}$`)
tab6_res <- tab6 |>
  select(Design, `$\\beta$`,
         `$\\lambda$`, `$\\kappa$`, `$\\eta$`,
         `$n_{\\max}$`, `$N_{H_0}$`, `$N_{H_1}$`, `$N_{H_0}^{P}$`, `$N_{H_1}^{P}$`,
         `$CP_{\\min}$`)
tn <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab6_n, k = c(0, 1, rep(0, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_6_n} Cumulative sample sizes of the Schlömer Brannath design (2) and the
		two proposed designs (4 with non-binding futility, 5 with binding futility)."
)
tb <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab6_b, k = c(0, 1, rep(2, 6)))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_6_b} Z-scale rejection boundaries of the three compared designs."
)
tr <- kbl(
  as_tibble(Map(OptimalGoldstandardDesigns:::fr, tab6_res, k = c(0, 1, 2, 1, 2, 0, 0, 0, 0, 0, 2))),
  format="latex",
  escape = FALSE,
  booktabs = TRUE,
  linesep = c("", "", "\\addlinespace"),
  align = 'c',
  caption = "\\label{tab_6_res} Operating characteristics and tuning parameters of the three compared designs."
)
cat(tn, tb, tr, sep = "\n\n")


### Figure 20 ###
# Conditional power plot
erglist <- list()
for (Z_TP in seq(-1.5, 2.5, by = 0.025)){
  for (Z_TC in seq(-1.5, 2.5, by = 0.025)){
    erglist[[length(erglist)+1]] <- list(Z_TP = Z_TP,
                                         Z_TC = Z_TC,
                                         cpower = OptimalGoldstandardDesigns:::calc_conditional_power(Z_TP,
                                                                                                      Z_TC,
                                                                                                      tab2[2,]$.design_object[[1]])
    )
  }
}

erg_df <- do.call(rbind.data.frame, erglist)
f_cpower <- ggplot(data = erg_df, aes(x=Z_TP, y=Z_TC, z=cpower)) +
  geom_contour_filled(breaks = c(-1, 0.025, seq(0.1, .9, 0.1), .99, 1.00001)) +
  labs(fill="Conditional power") +
  theme_pubclean() +
  scale_x_continuous(TeX("$Z_{1,TP}$"), expand = c(0,0.05), breaks = seq(-1, 5, 1)) +
  scale_y_continuous(TeX("$Z_{1,TC}$"), expand = c(0,0.05), breaks = seq(-1, 5, 1)) +
  theme(legend.position = "right",
        panel.grid.major.y = element_blank()) +
  scale_fill_viridis_d(alpha = 1, labels = c("[0%, 2.5%]",
                                             "(2.5%, 10%]",
                                             "(10%, 20%]",
                                             "(20%, 30%]",
                                             "(30%, 40%]",
                                             "(40%, 50%]",
                                             "(50%, 60%]",
                                             "(60%, 70%]",
                                             "(70%, 80%]",
                                             "(80%, 90%]",
                                             "(90%, 99%]",
                                             "early efficacy"),
                       option="D")
f_cpower
ggsave(filename = "figures/f_cpower.pdf",
       plot = f_cpower,
       device = "pdf",
       width = unit(7, "cm"), height = unit(5, "cm"))


