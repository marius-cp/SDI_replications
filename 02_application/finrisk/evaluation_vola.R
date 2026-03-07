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
devtools::install_github("marius-cp/SDI")
library(SDI)
source("../../00_functions/funs_plots.R")
set.seed(123)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# styling globals  ----
cols_models <-c("red", "blue", "purple", "coral")
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# evaluation mse ----
fcasts <- readRDS("data/emini_fcasts_vola.rds")

datse <-  fcasts %>% 
  pivot_longer(cols = HAR:GARCH, names_to = "model")
datse 

pdatse <- 
  datse %>% 
  group_by(model) %>% 
  summarise(
    MZdec(
      x = value, 
      y= rv5,
      S = function(x,y) 1*(x-y)^2
    )$dec
  )

# UNC SE 
pdatse$unc %>% unique() %>% round(2)

seq(min(pdatse$s),max(pdatse$s), length.out=8) %>% diff()
mcbdscplot_mse <- 
  ggplot()+
  geom_point(pdatse,mapping=aes(x=mcb, y=dsc, color=model), size=2)+
  coord_cartesian(ylim = c(8.5,11.5), xlim = c(0, .3)) +
  geom_abline(
    mapping=aes(
      intercept = isoline_intercept(
        s=sort(
          c(
            seq(min(pdatse$s),0, by= -0.3332847),
            seq(min(pdatse$s),max(pdatse$s), length.out=8),
            seq(max(pdatse$s),20, by= 0.3332847)
          ),
          decreasing = T),
        unc = na.omit(unique(pdatse$unc))[1]
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
        s=sort(c(min(pdatse$s),  max(pdatse$s)),decreasing = F),
        unc = na.omit(unique(pdatse$unc))[1]
      )$intercept, 
      slope = 1,
      label = isoline_intercept(
        s=sort(c(min(pdatse$s), max(pdatse$s)),decreasing = F),
        unc = na.omit(unique(pdatse$unc))[1]
      )$score %>% round(2)
    ),
    color = "darkgray",
    hjust = 0.9,
    size = 5,
    text_only = TRUE,
    boxcolour = NA,
    straight = TRUE
  )+
  geom_text(
    pdatse %>% mutate(dsc_adj = dsc+c(0,0,0,.25), mcb_adj=mcb+c(0,0,0,-.04)),
    mapping=aes(y=dsc_adj,x=mcb_adj),
    label= pdatse$model,
    size=5,
    hjust=.10, vjust=1.5
  )+
  scale_color_manual(values = cols_models)
mcbdscplot_mse

# Get unique models
ms <- datse$model %>% unique()

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
  # Get the volatility and RV_5min values for the current model
  X1 <- datse %>% filter(model == m) %>% pull(value)
  RV_5min <- datse %>% filter(model == m) %>% pull(rv5)
  
  # RELATIVE
  # Loop through rival models
  rivals <- ms[ms != m]
  for (r in rivals) {
    # Get the volatility for the rival model
    X2 <- datse %>% filter(model == r) %>% pull(value)
    
    # Calculate SDI between the current model and rival
    sdi <- SDI(
      X_1 = X1,      # RVar fcast
      X_2 = X2,      # RVar fcast
      Y = RV_5min,     # RVar 
      S = function(x,y) 1*(x-y)^2, 
      V = function(x,y) 1*2*(x-y),
      Spp = function(x, y) 1*2     
    )
    
    # Extract components from the SDI result and store in the data frames
    comparison_dataframes$s[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[1]
    comparison_dataframes$mcb[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[2]
    comparison_dataframes$dsc[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[3]
    comparison_dataframes$s_pval[m, r] <- sdi$asyvardm$pvals[1]
    comparison_dataframes$mcb_pval[m, r] <- sdi$asyvardm$pvals[2]
    comparison_dataframes$dsc_pval[m, r] <- sdi$asyvardm$pvals[3]
  }
  
  #ABSOLUT
  comparison_dataframes$mcb_null[m,m] <- sdi$mcbnulltest_X1$mcb_pval
  comparison_dataframes$MZ_classic[m,m] <-  NA #sdi$mcbnulltest_X1$mz_pval
  comparison_dataframes$dsc_null[m,m] <-  sdi$dscnulltest_X1$dsc_pval
}


facet_plot_mse <- plot_sdi_comparison_matrix(comparison_dataframes)

(free(mcbdscplot_mse, "label") | facet_plot_mse) +
  plot_layout(widths = c(1, 2))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# evaluation qlike ----
fcasts <- readRDS("data/emini_fcasts_vola.rds") 
fcasts

dat <-  fcasts %>% 
  pivot_longer(cols = HAR:GARCH, names_to = "model") 
dat

pdat <- 
  dat %>% 
  group_by(model) %>% 
  summarise(
    MZdec(
      x = value, 
      y= rv5,
      S = function(x, y) ((y / x) - log(y / x) - 1)*10 
    )$dec
  )

# UNC QLIKE
pdat$unc %>% unique() %>% round(3)

seq(min(pdat$s),max(pdat$s), length.out=6) %>% diff
mcbdscplot_qlike <- 
  ggplot()+
  geom_point(pdat,mapping=aes(x=mcb, y=dsc, color=model), size=2)+
  coord_cartesian(ylim = c(5.5, 7), xlim = c(0, 0.08)) +
  geom_abline(
    mapping=aes(
      intercept = isoline_intercept(
        s=sort(
          c(
            seq(min(pdat$s),0, by= -0.1753291),
            seq(min(pdat$s),max(pdat$s), length.out=6),
            seq(max(pdat$s),10, by= 0.1753291)
          ),
          decreasing = T),
        unc = na.omit(unique(pdat$unc))[1]
      )$intercept, 
      slope = 1
    ),
    color = "lightgray", alpha = 0.5,
    size = 0.5
  )+
  theme_classic()+
  themecompplot+
  xlab("MCB")+
  ylab("DSC")+
  geomtextpath::geom_labelabline(
    aes(
      intercept = isoline_intercept(
        s=sort(c(min(pdat$s),  max(pdat$s)),decreasing = F),
        unc = na.omit(unique(pdat$unc))[1]
      )$intercept, 
      slope = 1,
      label = isoline_intercept(
        s=sort(c(min(pdat$s), max(pdat$s)),decreasing = F),
        unc = na.omit(unique(pdat$unc))[1]
      )$score %>% round(2)
    ),
    color = "darkgray",
    hjust = 0.9,
    size = 5,
    text_only = TRUE,
    boxcolour = NA,
    straight = TRUE
  )+
  geom_text(
    pdat %>% mutate(dsc_adj = dsc+c(0,0,0,+0.1),mcb_adj = mcb+c(0,0,0,-0.009)),
    mapping=aes(y=dsc_adj,x=mcb_adj),
    label= pdat$model,
    size=5,
    hjust=0, vjust=1.5
  )+
  scale_color_manual(values = cols_models)

mcbdscplot_qlike

# Get unique models
ms <- dat$model %>% unique()

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
  # Get the volatility and RV_5min values for the current model
  X1 <- dat %>% filter(model == m) %>% pull(value)
  RV_5min <- dat %>% filter(model == m) %>% pull(rv5)
  
  # RELATIVE
  # Loop through rival models
  rivals <- ms[ms != m]
  for (r in rivals) {
    # Get the volatility for the rival model
    X2 <- dat %>% filter(model == r) %>% pull(value)
    
    # Calculate SDI between the current model and rival
    sdi <- SDI(
      X_1 = X1,      # RVar forecasts
      X_2 = X2,      # RVar forecasts
      Y = RV_5min,    # RVar forecasts
      S = function(x, y) ((y / x) - log(y / x) - 1)*10,  # QLIKE scoring function, scaled by 10
      V = function(x, y) ((x - y) / x^2)*10,             # First derivative of S with respect to x, scaled by 10
      Spp = function(x, y) (-(x - 2 * y) / x^3)*10       # Second derivative of S with respect to x, scaled by 10
    )
    
    # Extract components from the SDI result and store in the data frames
    comparison_dataframes$s[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[1]
    comparison_dataframes$mcb[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[2]
    comparison_dataframes$dsc[m, r] <- (sdi$asyvardm$dec_1 - sdi$asyvardm$dec_2)[3]
    comparison_dataframes$s_pval[m, r] <- sdi$asyvardm$pvals[1]
    comparison_dataframes$mcb_pval[m, r] <- sdi$asyvardm$pvals[2]
    comparison_dataframes$dsc_pval[m, r] <- sdi$asyvardm$pvals[3]
  }
  
  #ABSOLUT
  comparison_dataframes$mcb_null[m,m] <- sdi$mcbnulltest_X1$mcb_pval
  comparison_dataframes$MZ_classic[m,m] <- NA# sdi$mcbnulltest_X1$mz_pval
  comparison_dataframes$dsc_null[m,m] <-  sdi$dscnulltest_X1$dsc_pval
}

facet_plot_qlike <- plot_sdi_comparison_matrix(comparison_dataframes)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# table motivational example ----
# Extract model comparison results and display them
X1 <- dat %>% filter(model == "GARCH") %>% pull(value)
RV_5min <- dat %>% filter(model == "GARCH") %>% pull(rv5)
X2 <- dat %>% filter(model == "HAR") %>% pull(value)

# Calculate SDI between the current model and rival
sdi <- SDI(
  X_1 = X1,       
  X_2 = X2,       
  Y = RV_5min,      
  S = function(x, y) ((y / x) - log(y / x) - 1)*10,  # QLIKE scoring function, scaled by 10
  V = function(x, y) ((x - y) / x^2)*10,             # First derivative of S with respect to x, scaled by 10
  Spp = function(x, y) (-(x - 2 * y) / x^3)*10       # Second derivative of S with respect to x, scaled by 10
  
)

tibble(
  model = c("GARCH", "HAR"),
  MZpval = c(
    sdi$mcbnulltest_X1$mz_pval,
    sdi$mcbnulltest_X2$mz_pval
  )
) %>%
  # Combine additional decision metrics from the SDI analysis
  bind_cols(
    bind_rows(
      sdi$asyvardm$dec_1,
      sdi$asyvardm$dec_2
    )
  ) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  kbl(
    ., 
    format = "latex",
    booktabs = T, escape = T
  ) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# save ----
w<-12
h <- 11
row1 <- 
  (free(mcbdscplot_mse, "label") | facet_plot_mse) +
  plot_layout(widths = c(1, 2))
row1

row2 <- (free(mcbdscplot_qlike, "label") | facet_plot_qlike) +
  plot_layout(widths = c(1, 2))
row2


# Wrap each row so it's treated as a single element
wrapped_row1 <- wrap_elements(full = row1)
wrapped_row2 <- wrap_elements(full = row2)
wrapped_row1 <- wrap_elements(
  full = row1 +
    plot_annotation(title = "b) SE Score") &
    theme(plot.title = element_text(hjust = 0,size = 15))  # optional
)
wrapped_row2 <- wrap_elements(
  full = row2 +
    plot_annotation(title = "a) QLIKE Score") &
    theme(plot.title = element_text(hjust = 0,size = 15))
)
# 0.06232167 is the equal DSC pval for HAR/GARCH

# Stack and annotate just the rows
final_plot <- (wrapped_row2 / wrapped_row1) 
  plot_annotation(tag_levels = 'a', tag_suffix = ')')&
  theme(
    plot.tag.position = c(0.0125, .95),  # optional: fine-tune position (left, top)
    plot.tag = element_text(size = 15, margin = margin(b = 2)),  # Increase size here
  )

final_plot
w<-14
h <- 13
ggsave("plots/appl_vola_emini_comb.pdf", width = w, height = h,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/appl_vola_emini_comb.pdf",
  width = w, height = h,
  device = cairo_pdf
)
# sometimes problems by saving pdf, think it is due to phantom(0) for ggplot alignment in matrix
# running dev.off() and saving again helps