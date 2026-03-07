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
# globals
themeplot <- 
  theme(
    legend.position = "bottom",
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
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# plot data
fcastsvola <- readRDS("data/emini_fcasts_vola.rds") %>% 
  mutate(
    target ="Variance", level="Variance", 
    rv5=sqrt(rv5), HAR=sqrt(HAR),MIDAS=sqrt(MIDAS), GJR=sqrt(GJR), GARCH=sqrt(GARCH)
    ) #%>% select(-ret)
fcastsvar_1 <- readRDS("data/emini_fcasts_VaR_1.rds") %>% mutate(level="1% VaR", target ="VaR") #%>% select(-rv5)
fcastsvar_5 <- readRDS("data/emini_fcasts_VaR_5.rds") %>% mutate(level="5% VaR", target ="VaR") #%>% select(-rv5)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# plot
bind_rows(fcastsvola,fcastsvar_1, fcastsvar_5) %>% 
  filter(date >= "2019-01-01", date <= "2022-01-01") %>% 
  mutate(rel=ifelse(target=="VaR", ret, rv5)) %>% 
  select(date,rel,HAR,MIDAS,GJR,GARCH,HS,target,level) %>%
  pivot_longer(cols = HAR:HS, names_to = "model") %>% 
  mutate(
    # need to order the legend and facet enteries
    level = factor(level, levels = c("Variance", "1% VaR", "5% VaR")),
    target = factor(target, levels = c("Variance", "VaR")),
    
    model = factor(
      model,
      levels = c("HAR", "MIDAS","GARCH", "GJR","HS")  # HS last
    )
  ) %>% 
  filter(!(model=="HS" & target=="Variance")) %>% 
  ggplot(aes(x=date, y=value,color=level))+
  geom_line(aes(x=date, y=(rel)), color="darkgray")+
  geom_line(size=.75)+
  facet_grid2(
    model~target, 
    scales="free",
    independent = "y"
  )+
  scale_color_manual(
    values = c(
      "Variance"   = "blue",
      "1% VaR" = "coral",
      "5% VaR" = "purple"
    ),
    name = NULL
  )+
  theme_bw()+
  themeplot+
  xlab("Time")

w<-12
h <- 15
ggsave("plots/timeseries.pdf", width = w, height = h,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/timeseries.pdf", 
  width = w, height = h,  
  device = cairo_pdf
)

