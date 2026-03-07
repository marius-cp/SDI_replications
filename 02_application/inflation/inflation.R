rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(lubridate)
library(tidyverse)
library(kableExtra)
library(sandwich)
library(patchwork)
devtools::install_github("marius-cp/SDI")
source("../../00_functions/funs_plots.R")
library(SDI)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# load data ----
dat <- readRDS("./murphy_replication/inflation_mean.rds")

dat <- 
  readRDS("./murphy_replication/inflation_mean.rds") %>% 
    filter(
        stemp_rlz<="2020-12-01"
        )
dat$yq <- zoo::as.yearqtr(dat$stemp_rlz)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# TS plot ----
dat %>% 
  rename(
  "SPF Forecast  " = spf, 
  "Michigan Forecast  " = michigan, 
  "Inflation Rate" = rlz
  ) %>% 
  pivot_longer(cols = `SPF Forecast  `:`Inflation Rate`, values_to = "values", names_to = "id" ) %>% 
    mutate(id=factor(id, levels = c("SPF Forecast  ", "Michigan Forecast  ", "Inflation Rate"))) %>% 
    ggplot(aes(x=yq,y=values,color=id, linetype=id))+
    geom_line()+
    theme_bw()+
    theme(
        legend.position = "bottom",
        legend.title=element_blank(),
        legend.key.size = unit(1.5,"line"),
        legend.box.spacing = unit(0, "pt")
    )+
    #geom_vline(xintercept=as.Date("2022-08-01"))+
    xlab("")+
    ylab("Inflation Rate (percent)")+
    scale_color_manual(values=c("coral", "blue", "black"))+
    scale_linetype_manual(values=c("solid", "solid", 11))+
    zoo::scale_x_yearqtr(
        breaks = c(seq(zoo::as.yearqtr("1982Q4"), zoo::as.yearqtr("2020Q4"), length.out=8)),
        format = "%YQ%q"
        )
ggsave("plots/appl_inflation_timeseries.pdf", width = 7, height = 2.5,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/appl_inflation_timeseries.pdf", 
  width = 7, height = 2.5,  
  device = cairo_pdf
)

lm(dat$rlz ~dat$spf) %>% summary
lm(dat$rlz ~dat$michigan) %>% summary

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# components ---- 
infl <- 
    bind_rows(
        decomposition(
            x = dat$michigan,
            x_rc = fitted(lm(dat$rlz~dat$michigan)),
            y = dat$rlz, 
            S = function(x,y) (x-y)^2
            ) %>% mutate(id="Michigan"),
        
        decomposition(
            x = dat$spf,
            x_rc = fitted(lm(dat$rlz~dat$spf)),
            y = dat$rlz,
            S = function(x,y) (x-y)^2
        ) %>% mutate(id="SPF")
    )
infl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# p-values ----
vars <- asy_var_dm(
  X_1 = dat$michigan, 
  X_2 = dat$spf, 
  Y = dat$rlz, 
  S = function(x,y) (x-y)^2, 
  V = function(x,y) (x-y)*2
  )$asy_vars 
vars
pvals <- get_pval(dt = vars, tt=nrow(dat))

# UIIU
sdi <- SDI(
  X_1 = dat$michigan, 
  X_2 = dat$spf, 
  Y = dat$rlz, 
  S = function(x,y) (x-y)^2, 
  V = function(x,y) (x-y)*2, 
  Spp = function(x,y) 2
)
sdi$asyvardm$pvals$s_pval
sdi$mcbIUUI
sdi$dscIUUI

pvals %>% round(3) %>%
  pivot_longer(
    cols = s_pval:dsc_pval, 
    names_to = "test"
  ) %>% 
  mutate(order=c(3,1,2)) %>% 
  arrange(order) %>% 
  select(test, value) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  kbl(
    ., 
    format = "latex",
    booktabs = T, escape = T
  )# %>% 
  #readr::write_lines(file = "tabs/appl_inflation_pvalues.tex")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# plot components ---- 
comps <- 
ggplot()+
    theme_classic()+
    geom_point(data = infl, mapping = aes(y=dsc, x=mcb, color=id), size=2)+
    coord_cartesian(ylim = c(0, .8), xlim = c(0, .8)) +
    geom_abline(
        aes(
            intercept = isoline_intercept(
                s=sort(
                    c(
                        seq(0,1.462918, by=0.135145),
                        seq(1.462918,2.003499, length.out=5), 
                        seq(2.003499,5, by=0.135145)
                        ),
                    decreasing = T),
                #sort(c(seq(0,5,by=0.15)[-c(11,14)],c(2.003499, 1.462918)),decreasing = T),
                unc = na.omit(unique(infl$unc))[1]
            )$intercept, 
            slope = 1
        ),
        color = "lightgray", alpha = 0.5,
        size = 0.5
    )+
    geomtextpath::geom_labelabline(
        aes(
            intercept = isoline_intercept(
                s=sort(c(1.462918, 2.003499),decreasing = F),
                unc = na.omit(unique(infl$unc))[1]
            )$intercept, 
            slope = 1,
            label = isoline_intercept(
                s=sort(c(1.462918, 2.003499),decreasing = F),
                unc = na.omit(unique(infl$unc))[1]
            )$score %>% round(3)
            #intercept = intercept, slope = slope, label = Score
        ),
        color = "darkgray",
        hjust = 0.8,
        size = 5,#7*0.36,
        text_only = TRUE,
        boxcolour = NA,
        straight = TRUE
    )+
    geom_text(
        infl,
        mapping = aes(y=dsc,x=mcb),
        label= c("Michigan", "SPF"), #plotdat$id %>% unique()  ,
        size=5,
        hjust=1.2, vjust=.35)+
    xlab("MCB")+
    ylab("DSC")+
  theme(
  legend.position = "none",#c(0.95, 0.95),  # Position legend in upper left
  legend.title = element_blank(),
  legend.background = element_rect(fill = alpha('white', 0.5)),  # Optional: semi-transparent background
  legend.margin = margin(-5, 0, 0, 0),     # Remove margins around the legend
  legend.key.size = unit(6, "mm"),
  legend.text = element_text(
    colour = "darkgray",
    size = 15
  ),
  legend.spacing.x = unit(1, "mm"),
  axis.text.x = element_text(
    #angle = 45,
    #vjust = 1,
    #hjust=1,
    colour = "black",
    size = 15
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 15
  ),
  strip.text = element_text(
    size = 15#, face = "bold"
  ),
  axis.title=element_text(
    size=15
  ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
  #  strip.text.x = element_text(margin = margin(.5,.5,.5,.5, "mm")), 
  #   strip.text.y = element_text(margin = margin(.5,.5,.5,.5, "mm"))
)+
scale_color_manual(values=c("blue", "coral","black"))
comps
ggsave("plots/appl_inflation_componentplot.pdf", width = 3, height = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# MZ reg ----
mz_spf <- MZ.test.mean(haty = dat$spf, y = dat$rlz)
mz_spf_R2 <- summary(lm(dat$rlz ~ dat$spf))$r.squared
mz_michigan <- MZ.test.mean(haty = dat$michigan, y = dat$rlz)
mz_michigan_R2 <- summary(lm(dat$rlz ~ dat$michigan))$r.squared

mztab <- 
  rbind(
    c(
      "$\\beta_0$", 
      c(
        mz_spf$estimate[1],mz_michigan$estimate[1]
      ) %>% round(3)
    ),
    c(
      "$\\beta_1$", 
      c(
        mz_spf$estimate[2],mz_michigan$estimate[2]
      )  %>% round(3)
    ), 
    c(
      "$p$-value",
      c(
        mz_spf$Wald.pval,mz_michigan$Wald.pval
      ) %>% round(3)
    ), 
    c(
      "$R^2$",
      c(
        mz_spf_R2, mz_michigan_R2
      ) %>% round(3)
    )
  ) %>% 
  kbl(
    ., 
    format = "latex",
    col.names = c("",rep(c("SPF","Michigan"),1)), 
    booktabs = T, escape = T, digits = 3
  ) #%>% 
  #readr::write_lines(file = "tabs/appl_inflation_mztab.tex")
writeLines(mztab)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# matrix plot ----
dat
datse <-  dat %>% 
  #filter(Symbol==s) %>% 
  pivot_longer(cols = spf:michigan, names_to = "model") %>% 
  mutate(model=ifelse(model=="spf","SPF", "Michigan"))
datse

SDI(
  X_1 = dat$michigan, 
  X_2 = dat$spf, 
  Y = dat$rlz, 
  S = function(x,y) (x-y)^2, 
  V = function(x,y) (x-y)*2, 
  Spp = function(x,y) 2
)

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
  X1 <- datse %>% filter(model == m) %>% pull(value)
  rlz <- datse %>% filter(model == m) %>% pull(rlz)
  
  # RELATIVE
  # Loop through rival models
  rivals <- ms[ms != m]
  for (r in rivals) {
    X2 <- datse %>% filter(model == r) %>% pull(value)
    
    # Calculate SDI between the current model and rival
    sdi <- SDI(
      X_1 = X1,      
      X_2 = X2,      
      Y = rlz,   
      S = function(x,y) (x-y)^2, 
      V = function(x,y) 2*(x-y),
      Spp = function(x, y) 2      # Second derivative of S with respect to x
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

facet_plot_mse

comps | facet_plot_mse

w<-11
h <- 5.5
free(comps, "label")  | (
  (plot_spacer()/(facet_plot_mse))|plot_spacer()
  ) + 
  plot_layout(widths = c(1,.5))
ggsave("plots/appl_infl.pdf", width = w, height = h,  device = cairo_pdf)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/appl_infl.pdf", 
  width = w, height = h,  
  device = cairo_pdf
)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

