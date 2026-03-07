rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(lubridate)
library(tidyverse)
library(sandwich)
library(patchwork)
library(kableExtra)
library(rugarch)
library(doParallel)
library(ggh4x)
library(gtable)
library(ggtext)
library(gridExtra)
library(grid)
library(knitr)
devtools::install_github("marius-cp/SDI")
library(SDI)
source("../../00_functions/funs_plots.R")
source("../../00_functions/funs_backtest.R")
set.seed(123)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# styling globals ----
cols_models <-c("red", "blue", "purple", "black", "coral")
themecompplot <- 
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha('white', 0.5)),  
    legend.margin = margin(-5, 0, 0, 0),     
    legend.key.size = unit(6, "mm"),
    legend.text = element_text(
      colour = "black",
      size = 15
    ),
    legend.spacing.x = unit(1, "mm"),
    axis.text.x = element_text(colour = "black",size = 15),
    axis.text.y = element_text(colour = "black",size = 15),
    strip.text = element_text(colour = "black",size = 15),
    axis.title= element_text(colour = "black",size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
# themecompplot <- 
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = alpha('white', 0.5)),  
#     legend.margin = margin(-5, 0, 0, 0),     
#     legend.key.size = unit(6, "mm"),
#     legend.text = element_text(
#       colour = "black",
#       size = 15
#     ),
#     legend.spacing.x = unit(1, "mm"),
#     axis.text.x = element_text(colour = "black",size = 15),
#     axis.text.y = element_text(colour = "black",size = 15),
#     strip.text = element_text(colour = "black",size = 15),
#     axis.title= element_text(colour = "black",size = 15),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# scores ----
scoreingfunction = function(x,y, alpha) (1 * (y < x) - alpha) * (x - y)
identificationfunction = function(x,y, alpha) (1 * (y < x) - alpha)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1% VaR ----
alpha = 0.01
fcasts <- readRDS("data/emini_fcasts_VaR_1.rds")

datpinball1 <-  fcasts %>% pivot_longer(cols = HAR:HS, names_to = "model")

pdatoinball1 <- 
  datpinball1 %>% 
  group_by(model) %>% 
  summarise(
    qMZdec(
      x = value, 
      y= ret,
      S = scoreingfunction, 
      alpha=alpha
    )$dec
  ) %>% 
  mutate(across(where(is.numeric), ~ 100 * .x))

seq(min(pdatoinball1$s),max(pdatoinball1$s), length.out=8) %>% diff()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## component plot ----
mcbdscplot1 <- 
  ggplot()+
  geom_point(pdatoinball1,mapping=aes(x=mcb, y=dsc, color=model), size=2)+
  coord_cartesian(ylim = c(0.03,2), xlim = c(0, .4)) +
  geom_abline(
    mapping=aes(
      intercept = isoline_intercept(
        s=sort(
          c(
            seq(min(pdatoinball1$s),0, by= -0.1702278),
            seq(min(pdatoinball1$s),max(pdatoinball1$s), length.out=8),
            seq(max(pdatoinball1$s),20, by= 0.1702278)
          ),
          decreasing = T),
        unc = na.omit(unique(pdatoinball1$unc))[1]
      )$intercept, 
      slope = 1
    ),
    color = "lightgray", alpha = 0.5,
    size = 0.5
  )+
  theme_classic()+
  themecompplot +
  xlab("MCB")+
  ylab("DSC")+
  geomtextpath::geom_labelabline(
    aes(
      intercept = isoline_intercept(
        s=sort(c(min(pdatoinball1$s),  max(pdatoinball1$s)),decreasing = F),
        unc = na.omit(unique(pdatoinball1$unc))[1]
      )$intercept, 
      slope = 1,
      label = isoline_intercept(
        s=sort(c(min(pdatoinball1$s), max(pdatoinball1$s)),decreasing = F),
        unc = na.omit(unique(pdatoinball1$unc))[1]
      )$score %>% round(2)
      #intercept = intercept, slope = slope, label = Score
    ),
    color = "darkgray",
    hjust = 0.9,
    size = 5,
    text_only = TRUE,
    boxcolour = NA,
    straight = TRUE
  )+
  geom_text(
    #pdatoinball1 %>% mutate(dsc_adj = dsc+c(0,0,.4,0,0,0)),
    pdatoinball1 %>% mutate(dsc_adj = dsc+c(0,0,0,0,0.15), mcb_adj=mcb+c(0,0,0,0,-0.025)),
    mapping=aes(y=dsc_adj,x=mcb_adj),
    label= pdatoinball1$model,
    size=5,
    hjust=.10, vjust=1.5
  )+
  scale_color_manual(values = cols_models)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## matrix plot ----

# Get unique models
ms <- datpinball1$model %>% unique()

# Initialize data frames for s, mcb, dsc, and their p-values
comparison_dataframes <- list(
  s = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  mcb = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  dsc = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  s_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  mcb_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  dsc_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  mcb_null = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  MZ_classic = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  dsc_null = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms)))
)

# Loop through each model
for (m in ms) {
  X1 <- datpinball1 %>% filter(model == m) %>% pull(value)
  ret <- datpinball1 %>% filter(model == m) %>% pull(ret)
  
  # RELATIVE
  # Loop through rival models
  rivals <- ms[ms != m]
  for (r in rivals) {
    # Get the volatility for the rival model
    X2 <- datpinball1 %>% filter(model == r) %>% pull(value)
    
    # Calculate SDI between the current model and rival
    sdi <- qSDI(
      X_1 = X1,      #  fcast
      X_2 = X2,      #  fcast
      Y = ret,     #
      S = scoreingfunction,
      V = identificationfunction,
      alpha=alpha
    )
    
    # Extract components from the SDI result and store in the data frames
    comparison_dataframes$s[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[1] *100
    comparison_dataframes$mcb[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[2] *100
    comparison_dataframes$dsc[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[3] *100
    comparison_dataframes$s_pval[m, r] <- sdi$asyvardm$pvals[1]
    comparison_dataframes$mcb_pval[m, r] <- sdi$asyvardm$pvals[2]
    comparison_dataframes$dsc_pval[m, r] <- sdi$asyvardm$pvals[3]
  }
  
  #ABSOLUT
  comparison_dataframes$mcb_null[m,m] <- sdi$mcbnulltest_X1$Wald.pval.MCB
  comparison_dataframes$MZ_classic[m,m] <-  NA#sdi$mcbnulltest_X1$mz_pval
  comparison_dataframes$dsc_null[m,m] <-  sdi$dscnulltest_X1$Wald.pval.DSC
}


facet_plot1 <- plot_sdi_comparison_matrix_q(comparison_dataframes)


row1 <- (free(mcbdscplot1, "label") | facet_plot1) + plot_layout(widths = c(1, 2))
row1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## backtest table ----

backtesttab1 <- 
  datpinball1 %>%
  filter(model %in% c("HAR", "MIDAS", "GJR", "GARCH", "HS")) %>%
  mutate(alpha=alpha,hit_demeaned = as.numeric(ret < value) - alpha) %>%
  arrange(model, date) %>%
  group_by(model) %>%
  mutate(
    hit_lag1 = lag(hit_demeaned, 1),
    hit_lag2 = lag(hit_demeaned, 2),
    hit_lag3 = lag(hit_demeaned, 3),
    hit_lag4 = lag(hit_demeaned, 4)
  ) %>%
  summarise(
    backtest = list(backtesttab_wrapper(cur_data())),
    .groups = "drop"
  ) %>%
  tidyr::unnest(backtest) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  data.frame()

backtesttab1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## mot ex tab ----
bind_cols(
  pdatoinball1 %>%  filter(model %in% c("GARCH", "HAR", "HS")),
  data.frame(
    MZp =c(
      comparison_dataframes$mcb_null[4,4],#garch
      comparison_dataframes$mcb_null[1,1],#har
      comparison_dataframes$mcb_null[5,5]#hs
    ) %>% round(3)
  )
  ) %>% 
  arrange(MZp) %>% 
  select(model, MZp, s, mcb, dsc, unc) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  data.frame()

# pval  DM test
comparison_dataframes$s_pval
comparison_dataframes$s_pval[4,1] #HAR vs GARCH 0.486
comparison_dataframes$s_pval[5,1] #HAR vs HS 0.000
comparison_dataframes$s_pval[4,5] #GARCH vs HS 0.0004

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 5% VaR ----
alpha =0.05
fcasts <- readRDS("data/emini_fcasts_VaR_5.rds")

datpinball5 <-  fcasts %>% pivot_longer(cols = HAR:HS, names_to = "model")

pdatoinball5 <- 
  datpinball5 %>% 
  group_by(model) %>% 
  summarise(
    qMZdec(
      x = value, 
      y= ret,
      S = scoreingfunction, 
      alpha=alpha
    )$dec
  ) %>% 
  mutate(across(where(is.numeric), ~ 100 * .x))

seq(min(pdatoinball5$s),max(pdatoinball5$s), length.out=8) %>% diff()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## component plot ----
mcbdscplot5 <- 
  ggplot()+
  geom_point(pdatoinball5,mapping=aes(x=mcb, y=dsc, color=model), size=2)+
  coord_cartesian(ylim = c(1,4), xlim = c(0, .21)) +
  geom_abline(
    mapping=aes(
      intercept = isoline_intercept(
        s=sort(
          c(
            seq(min(pdatoinball5$s),0, by= -0.3897666),
            seq(min(pdatoinball5$s),max(pdatoinball5$s), length.out=8),
            seq(max(pdatoinball5$s),20, by= 0.3897666)
          ),
          decreasing = T),
        unc = na.omit(unique(pdatoinball5$unc))[1]
      )$intercept, 
      slope = 1
    ),
    color = "lightgray", alpha = 0.5,
    size = 0.5
  )+
  theme_classic()+
  themecompplot +
  xlab("MCB")+
  ylab("DSC")+
  geomtextpath::geom_labelabline(
    aes(
      intercept = isoline_intercept(
        s=sort(c(min(pdatoinball5$s),  max(pdatoinball5$s)),decreasing = F),
        unc = na.omit(unique(pdatoinball5$unc))[1]
      )$intercept, 
      slope = 1,
      label = isoline_intercept(
        s=sort(c(min(pdatoinball5$s), max(pdatoinball5$s)),decreasing = F),
        unc = na.omit(unique(pdatoinball5$unc))[1]
      )$score %>% round(2)
      #intercept = intercept, slope = slope, label = Score
    ),
    color = "darkgray",
    hjust = 0.9,
    size = 5,
    text_only = TRUE,
    boxcolour = NA,
    straight = TRUE
  )+
  geom_text(
    #pdatoinball5 %>% mutate(dsc_adj = dsc+c(0,0,.4,0,0,0)),
    pdatoinball5 %>% mutate(dsc_adj = dsc+c(0,0,0,0,0.2), mcb_adj=mcb+c(0,0,0,0,0)),
    mapping=aes(y=dsc_adj,x=mcb_adj),
    label= pdatoinball5$model,
    size=5,
    hjust=.10, vjust=1.5
  )+
  scale_color_manual(values = cols_models)
mcbdscplot5

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## matrix plot ----

# Get unique models
ms <- datpinball5$model %>% unique()

# Initialize data frames for s, mcb, dsc, and their p-values
comparison_dataframes <- list(
  s = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  mcb = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  dsc = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  s_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  mcb_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))),
  dsc_pval = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  mcb_null = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  MZ_classic = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms))), 
  dsc_null = as.data.frame(matrix(NA, nrow = length(ms), ncol = length(ms), dimnames = list(ms, ms)))
)

# Loop through each model
for (m in ms) {
  X1 <- datpinball5 %>% filter(model == m) %>% pull(value)
  ret <- datpinball5 %>% filter(model == m) %>% pull(ret)
  
  # RELATIVE
  # Loop through rival models
  rivals <- ms[ms != m]
  for (r in rivals) {
    # Get the volatility for the rival model
    X2 <- datpinball5 %>% filter(model == r) %>% pull(value)
    
    # Calculate SDI between the current model and rival
    sdi <- qSDI(
      X_1 = X1,      #  fcast
      X_2 = X2,      #  fcast
      Y = ret,     #
      S = scoreingfunction,
      V = identificationfunction,
      alpha=alpha
    )
    
    # Extract components from the SDI result and store in the data frames
    comparison_dataframes$s[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[1] *100
    comparison_dataframes$mcb[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[2] *100
    comparison_dataframes$dsc[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[3] *100
    comparison_dataframes$s_pval[m, r] <- sdi$asyvardm$pvals[1]
    comparison_dataframes$mcb_pval[m, r] <- sdi$asyvardm$pvals[2]
    comparison_dataframes$dsc_pval[m, r] <- sdi$asyvardm$pvals[3]
  }
  
  #ABSOLUT
  comparison_dataframes$mcb_null[m,m] <- sdi$mcbnulltest_X1$Wald.pval.MCB
  comparison_dataframes$MZ_classic[m,m] <-  NA#sdi$mcbnulltest_X1$mz_pval
  comparison_dataframes$dsc_null[m,m] <-  sdi$dscnulltest_X1$Wald.pval.DSC
}


facet_plot5 <- plot_sdi_comparison_matrix_q(comparison_dataframes)


row2 <- (free(mcbdscplot5, "label") | facet_plot5) + plot_layout(widths = c(1, 2))
row2

## backtest table ----
backtesttab5 <- 
  datpinball5 %>%
  filter(model %in% c("HAR", "MIDAS", "GJR", "GARCH", "HS")) %>%
  mutate(alpha=alpha,hit_demeaned = as.numeric(ret < value) - alpha) %>%
  arrange(model, date) %>%
  group_by(model) %>%
  mutate(
    hit_lag1 = lag(hit_demeaned, 1),
    hit_lag2 = lag(hit_demeaned, 2),
    hit_lag3 = lag(hit_demeaned, 3),
    hit_lag4 = lag(hit_demeaned, 4)
  ) %>%
  summarise(
    backtest = list(backtesttab_wrapper(cur_data())),
    .groups = "drop"
  ) %>%
  tidyr::unnest(backtest) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  data.frame()

backtesttab5

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# final plot ----
wrapped_row1 <- wrap_elements(
  full = row1 +
    plot_annotation(title = "a) VaR 1%") &
    theme(plot.title = element_text(hjust = 0,size = 15))  # optional
)
wrapped_row2 <- wrap_elements(
  full = row2 +
    plot_annotation(title = "b) VaR 5%") &
    theme(plot.title = element_text(hjust = 0,size = 15))
)

final_plot <- wrapped_row1 / wrapped_row2
final_plot
w<-14
h <- 13
ggsave("plots/appl_VaR_emini_comb.pdf", width = w, height = h,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/appl_VaR_emini_comb.pdf",
  width = w, height = h,
  device = cairo_pdf
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# final table ----
tab_df <-
  bind_rows(
    backtesttab1 %>% mutate(ForecastType = "$1\\%$-Quantile"),
    backtesttab5 %>% mutate(ForecastType = "$5\\%$-Quantile")
  ) %>%
  select(ForecastType, model, UC, Basel, CC,NZ,  VQR, DQ,DQv, Baselc, relHits) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% 
  mutate(
    Basel = cell_spec(
      Basel,
      format = "latex",
      color  = "black",
      bold   = FALSE,
      background = case_when(
        Baselc == "green"  ~ "ForestGreen!15",
        Baselc == "yellow" ~ "Goldenrod!20",
        Baselc == "red"    ~ "Red!15",
        TRUE               ~ "gray!10"
      )
    )
  ) %>%
  select(-Baselc)
# 1) enforce model order + sort within each ForecastType block
tab_df <- tab_df %>%
  mutate(
    model = factor(model, levels = c("HS", "GARCH", "GJR", "HAR", "MIDAS"))
  ) %>%
  arrange(ForecastType, model)

# 2) create the multirow labels robustly 
tab_df <- tab_df %>%
  group_by(ForecastType) %>%
  mutate(
    ForecastType = if_else(
      row_number() == 1,
      paste0("\\multirow{", n(), "}{*}{", ForecastType, "}"),
      ""
    )
  ) %>%
  ungroup()

colnames(tab_df) <- c("", "Model $i$", "UC", "Basel", "CC", "NZ" , "VQR", "DQ", "DQX", "relative")

tab_latex <-
  tab_df %>%
  kbl(
    format="latex", 
    booktabs=TRUE, 
    escape=FALSE, 
    digits=3,
    caption="Backtest results."
    ) %>%
  add_header_above(
    c(
      " " = 2, 
      "p-values backtests" = 7, 
      "hit frequencies" = 1
      ),
    escape = FALSE
    ) %>%
  row_spec(
    5, 
    extra_latex_after="\\addlinespace\\hline\\addlinespace\n"
    )

tab_latex

