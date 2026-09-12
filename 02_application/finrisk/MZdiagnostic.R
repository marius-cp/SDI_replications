rm(list = ls())
if (!interactive()) grDevices::pdf(NULL)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
rm(script_arg)
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
library(SDI)
source("../../00_functions/funs_plots.R")
set.seed(123)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# globals
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
# Vola  ----
# Read forecasts dataset (wide format: one column per model HAR...GARCH)
fcasts <- readRDS("data/emini_fcasts_vola.rds") 
fcasts

# Build per-model recalibration and residuals
pdatse <- 
  fcasts %>% 
  pivot_longer(cols = HAR:GARCH, names_to = "model") %>% 
  mutate(
    model = factor(model, levels = c("HAR", "MIDAS", "GARCH", "GJR"))
  ) %>% 
  group_by(model) %>% 
  summarise(
    # NOTE: summarise() is *usually* used to collapse to 1 row per group.
    # Here you're returning full-length vectors (date/value/rv5/ret),
    # which technically works, but is unconventional. mutate() is the
    # more standard choice for this.
    
    date  = date,     # carry date through (vector)
    value = value,    # model forecast (vector)
    rv5   = rv5,      # realized volatility target (vector)
    ret   = ret,      # returns (vector)
    
    # MZdec recalibration:
    # - x = forecast, y = realized value
    # - S is squared error loss
    # - you take the recalibrated x (x_rc), i.e. the Mincer–Zarnowitz-type fit
    MZfit = MZdec(
      x = value, 
      y = rv5,
      S = function(x,y) 1*(x-y)^2
    )$x_rc
  ) %>% 
  # residual = realized - fitted (so positive means underprediction)
  mutate(
    res = (rv5 - MZfit), 
    model = as.factor(model)
    )

seplot <- 
ggplot(pdatse, aes(x = value, y = res)) +
  geom_point(alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  scale_x_log10()+
  coord_cartesian(ylim=c(-5,5))+
  theme_bw()+
  themediagplot+
  ylab("Residual")+
  xlab("Forecast")
seplot

# Build per-model recalibration and residuals
pdatqlike <- 
  fcasts %>% 
  pivot_longer(cols = HAR:GARCH, names_to = "model") %>% 
  mutate(
    model = factor(model, levels = c("HAR", "MIDAS", "GARCH", "GJR"))
  ) %>% 
  group_by(model) %>% 
  summarise(
    # NOTE: summarise() is *usually* used to collapse to 1 row per group.
    # Here you're returning full-length vectors (date/value/rv5/ret),
    # which technically works, but is unconventional. mutate() is the
    # more standard choice for this.
    
    date  = date,     # carry date through (vector)
    value = value,    # model forecast (vector)
    rv5   = rv5,      # realized volatility target (vector)
    ret   = ret,      # returns (vector)
    
    MZfit = MZdec(
      x = value, 
      y = rv5,
      S = function(x, y) ((y / x) - log(y / x) - 1)*10
    )$x_rc
  ) %>% 
  # residual = realized - fitted (so positive means underprediction)
  mutate(
    res = (rv5 - MZfit), 
    model = as.factor(model)
  )

qlikeplot <- 
ggplot(pdatqlike, aes(x = value, y = res)) +
  geom_point(alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  scale_x_log10()+
  coord_cartesian(ylim=c(-5,5))+
  theme_bw()+
  ylab("residual")+
  xlab("forecast")+
  theme_bw()+
  themediagplot+
  ylab("Residual")+
  xlab("Forecast")

wrapped_row1 <- wrap_elements(
  full = seplot +
    plot_annotation(title = "b) Mincer-Zarnowitz Regression with SE Score") &
    theme(plot.title = element_text(hjust = 0,size = 15))  # optional
)
wrapped_row2 <- wrap_elements(
  full = qlikeplot +
    plot_annotation(title = "a) Mincer-Zarnowitz Regression with QLIKE Score") &
    theme(plot.title = element_text(hjust = 0,size = 15))
)

# Stack and annotate just the rows
final_plot <- (wrapped_row2 / wrapped_row1) 
final_plot

w <- 12
h <- 8
ggsave("plots/MZdiagnostics_vola.pdf", width = w, height = h,  device = cairo_pdf)

# Deliberately disabled: paper-folder copies must be made manually.
# ggsave(
#   "/path/to/Overleaf/fig/MZdiagnostics_vola.pdf",
#   width = w, height = h,
#   device = cairo_pdf
# )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# VaR ----
## 1% ----
fcasts1 <- readRDS("data/emini_fcasts_VaR_1.rds") 
fcasts1

scoreingfunction = function(x,y, alpha) (1 * (y < x) - alpha) * (x - y)
identificationfunction = function(x,y, alpha) (1 * (y < x) - alpha)

# Build per-model recalibration and residuals
pdat1 <- 
  fcasts1 %>% 
  pivot_longer(cols = HAR:HS, names_to = "model") %>% 
  mutate(
    model = factor(model, levels = c("HAR", "MIDAS", "GARCH", "GJR", "HS"))
  ) %>% 
  group_by(model) %>% 
  summarise(
    date  = date,     # carry date through (vector)
    value = value,    # model forecast (vector)
    ret   = ret,      # returns (vector)
    MZfit = qMZdec(
      x = value, 
      y = ret,
      S = scoreingfunction, alpha = 0.01
    )$x_rc
  ) %>% 
  # residual = realized - fitted (so positive means underprediction)
  mutate(
    genres = identificationfunction(y=ret, x=MZfit, alpha = 0.01),
    model = as.factor(model)
  )

p1 <- 
ggplot(pdat1, aes(x = value, y = genres)) +
  geom_point(alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  #scale_x_log10()+
  coord_cartesian(ylim=c(-.02,1))+
  theme_bw()+
  ylab("residual")+
  xlab("forecast")+
  theme_bw()+
  themediagplot+
  ylab("Generalized Residual") +
  xlab("Forecast")
p1

# y-axis break
shift <- .8
p1 <-
  ggplot(
    pdat1 %>% mutate(genresmod = ifelse(genres==0.99,0.99-shift, genres)),
    aes(x = value, y = genres)) +
  geom_point(mapping=aes(x = value, y = genresmod),alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  #scale_x_log10()+
  #coord_cartesian(ylim=c(-.125,1-shift))+
  theme_bw()+
  ylab("residual")+
  xlab("forecast")+
  theme_bw()+
  themediagplot+
  ylab("Generalized Residual") +
  xlab("Forecast")+
  scale_y_continuous(
    breaks = c(-.1, -0.05, 0, 0.05, .1, .15,.2),
    labels = c("-0.1", "-0.05", "0", "0.05",
               "<span style='color:red;'>//</span>",
               "0.9", "1")
  ) +
  theme(
    axis.text.y = ggtext::element_markdown()
  )+
  geom_abline(intercept = .1, slope = 0, color ="red", linetype="dotted")

p1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## 5% ----
fcasts5 <- readRDS("data/emini_fcasts_VaR_5.rds") 
fcasts5

scoreingfunction = function(x,y, alpha) (1 * (y < x) - alpha) * (x - y)
identificationfunction = function(x,y, alpha) (1 * (y < x) - alpha)

# Build per-model recalibration and residuals
pdat5 <- 
  fcasts5 %>% 
  pivot_longer(cols = HAR:HS, names_to = "model") %>% 
  mutate(
    model = factor(model, levels = c("HAR", "MIDAS", "GARCH", "GJR", "HS"))
  ) %>% 
  group_by(model) %>% 
  summarise(
    date  = date,     # carry date through (vector)
    value = value,    # model forecast (vector)
    ret   = ret,      # returns (vector)
    MZfit = qMZdec(
      x = value, 
      y = ret,
      S = scoreingfunction, alpha = 0.05
    )$x_rc
  ) %>% 
  # residual = realized - fitted (so positive means underprediction)
  mutate(
    genres = identificationfunction(y=ret, x=MZfit, alpha = 0.05),
    model = as.factor(model)
  )

p5 <- 
  ggplot(pdat5, aes(x = value, y = genres)) +
  geom_point(alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  #scale_x_log10()+
  coord_cartesian(ylim=c(-.075,1))+
  theme_bw()+
  ylab("residual")+
  xlab("forecast")+
  theme_bw()+
  themediagplot+
  ylab("Generalized Residual") +
  xlab("Forecast")
p5

# y-axis break
shift <- .5
p5 <-
ggplot(
  pdat5 %>% mutate(genresmod = ifelse(genres==0.95,0.95-shift, genres)),
  aes(x = value, y = genres)) +
  geom_point(mapping=aes(x = value, y = genresmod),alpha = 0.25, color="gray") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ model, scales = "free_x") +
  geom_smooth(method="loess", fill="blue", alpha=.25)+
  #scale_x_log10()+
  coord_cartesian(ylim=c(-.125,1-shift))+
  theme_bw()+
  ylab("residual")+
  xlab("forecast")+
  theme_bw()+
  themediagplot+
  ylab("Generalized Residual") +
  xlab("Forecast")+
  scale_y_continuous(
    breaks = c(-.1, 0, .1, .2, .3, .4, .5),
    labels = c("-0.1", "0", "0.1", "0.2",
               "<span style='color:red;'>//</span>",
               "0.9", "1")
  ) +
  theme(
    axis.text.y = ggtext::element_markdown()
  )+
  geom_abline(intercept = .3, slope = 0, color ="red", linetype="dotted")

wrapped_row1 <- wrap_elements(
  full = p1 +
    plot_annotation(title = "a) Quantile Mincer-Zarnowitz Regression for 1% VaR") &
    theme(plot.title = element_text(hjust = 0,size = 15))  # optional
)
wrapped_row2 <- wrap_elements(
  full = p5 +
    plot_annotation(title = "b) Quantile Mincer-Zarnowitz Regression for 5% VaR") &
    theme(plot.title = element_text(hjust = 0,size = 15))  # optional
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Stack and annotate just the rows
final_plot <- (wrapped_row1 / wrapped_row2) 
final_plot
w <- 14
h <- 8
ggsave("plots/MZdiagnostics_VaR.pdf", width = w, height = h,  device = cairo_pdf)

# Deliberately disabled: paper-folder copies must be made manually.
# ggsave(
#   "/path/to/Overleaf/fig/MZdiagnostics_VaR.pdf",
#   width = w, height = h,
#   device = cairo_pdf
# )
