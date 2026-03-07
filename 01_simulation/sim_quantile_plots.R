rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")
lsize <- 1.05

# Main Text Simulation Figure -----
# dat <- readRDS("data/sim_q_parameterized.rds")
files <- paste0("data/sim_q_parameterized_part", 1:4, ".rds")
dat <- bind_rows(lapply(files, readRDS))

# Rejection rates  MCB test
rrmcb <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(1,2),
    type%in%c(
      #"pval", 
      #"bonf_pval", 
      "IUUI"
      #"trueX1", 
      #"trueX2"
    )
  ) %>% 
  dplyr::select(-s,-unc, -dsc, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k, alpha) %>% 
  summarise(
    rr = mean(mcb<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    k=unique(k), # changed from j to k
    alpha=unique(alpha)
  ) 

# Rejection rates DM test MCB simulation
rrs_mcb <- 
  dat %>% 
  filter(
    setup %in% c(1,2),
    type%in%c("pval")
  ) %>% 
  dplyr::select(-unc, -dsc,-mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k, alpha) %>% 
  summarise(
    rr = mean(s<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j), 
    alpha=unique(alpha)
  ) %>% 
  mutate(
    type="DM",
    panel="mcb"
  )

# Rejection rates  MCB test
rrdsc <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(3,4),
    type%in%c(
      "IUUI"#, 
      #"trueX1", "trueX2"
    )
  ) %>% 
  dplyr::select(-s,-unc, -mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k, alpha) %>% 
  summarise(
    rr = mean(dsc<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j), 
    alpha=unique(alpha)
  ) 


# Rejection rates DM test DSC simulation
rrs_dsc <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(3,4),
    type%in%c("pval")
  ) %>% 
  dplyr::select(-unc, -dsc,-mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k, alpha) %>% 
  summarise(
    rr = mean(s<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    k=unique(k), # changed from j to k
    alpha=unique(alpha)
  ) %>% 
  mutate(
    type="DM",
    panel="dsc"
  )



datall <- 
  bind_rows(
    rrmcb  %>% mutate(panel = "mcb", type = "proposed test on MCB difference"),
    rrs_mcb %>% mutate(panel = "mcb", type = "DM"),
    rrdsc  %>% mutate(panel = "dsc", type = "proposed test on DSC difference"),
    rrs_dsc %>% mutate(panel = "dsc", type = "DM")
  )%>% 
  mutate(
    type=factor(type),
    type =recode_factor(
      type,
      "proposed test on MCB difference" = "proposed test on MCB difference        ",
      "proposed test on DSC difference" = "proposed test on DSC difference        ",
      "DM"                     = "DM test        ",
      .ordered = TRUE   # this fixes the legend order to input order
    )
  )  %>%
  mutate(alpha_lab = paste0("alpha == ", alpha))

colors <- c("blue", "coral", alpha("black",0.5))
lines <- c("solid", "solid", "dotted")
datall %>% 
  ggplot(aes(y = rr, x = k, color = type,linetype = type)) +
  geom_abline(slope = 0, intercept = .1, color = "gray") +
  geom_line(size=lsize) +
  facet_grid(alpha_lab ~ mcbinfo + dscinfo, labeller = label_parsed) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(25, "mm"),
    legend.key.size = unit(25, "mm"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  xlab(expression(italic(k))) +
  ylab("empirical rejection rate")+
  scale_linetype_manual(
    values=lines
  )+
  scale_color_manual(    
    values = colors
  )
ggsave("plots/sim_q.pdf", height=10, width=10)
ggsave("/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/sim_q.pdf", height=10, width=10)

