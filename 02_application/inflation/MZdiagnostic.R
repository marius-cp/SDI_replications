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
themediagplot <- 
  theme(
    legend.position = "none",
    legend.text = element_text(
      colour = "black",
      size = 12
    ),
    legend.spacing.x = unit(1, "mm"),
    axis.text.x = element_text(colour = "black",size = 12),
    axis.text.y = element_text(colour = "black",size = 12),
    strip.text = element_text(colour = "black",size = 12),
    axis.title= element_text(colour = "black",size = 12)
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
fcasts <- readRDS("./murphy_replication/inflation_mean.rds") %>% filter(stemp_rlz<="2020-12-01")
fcasts

# Build per-model recalibration and residuals
pdatse <- 
  fcasts %>% 
  rename(
    "SPF Forecast  " = spf, 
    "Michigan Forecast  " = michigan, 
  ) %>% 
  pivot_longer(cols = "SPF Forecast  ":"Michigan Forecast  ", names_to = "model") %>% 
  group_by(model) %>% 
  summarise(
    dt = dt,     # carry date through (vector)
    value = value,    # model forecast (vector)
    rlz   = rlz,      
    MZfit = MZdec(
      x = value, 
      y = rlz,
      S = function(x,y) 1*(x-y)^2
    )$x_rc
  ) %>% 
  mutate(
    res = (rlz - MZfit), 
    model = as.factor(model)
    )


seplot <- 
ggplot(pdatse, aes(x = value, y = res)) +
  geom_point(alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  #scale_x_log10()+
  coord_cartesian(ylim=c(-5,5))+
  theme_bw()+
  themediagplot+
  ylab("Residual")+
  xlab("Forecast")

seplot

w <- 8
h <- 4
ggsave("plots/MZdiagnostics_inflation.pdf", width = w, height = h,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/MZdiagnostics_inflation.pdf",
  width = w, height = h,
  device = cairo_pdf
)



