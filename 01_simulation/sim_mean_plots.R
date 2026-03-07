rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(ggh4x)
library(tidyr)
library(ggnewscale)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")

K <- seq(0.0,.5, length.out=11)
dat <- readRDS("data/sim_m_parameterized.rds")


# Rejection rates  MCB test
rrmcb <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(1:6),
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
  group_by(setup,tt,type, k) %>% 
  summarise(
    rr = mean(mcb<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j), 
    k=unique(k)
  ) 
rrmcb

# Rejection rates DM test MCB simulation
rrs_mcb <- 
  dat %>% 
  filter(
    setup %in% c(
      2,4,6
    ),
    type%in%c("pval")
  ) %>% 
  dplyr::select(-unc, -dsc,-mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k) %>% 
  summarise(
    rr = mean(s<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j),     
    k=unique(k)
  ) %>% 
  mutate(
    type="DM",
    panel="mcb"
  )

# Rejection rates  MCB test
rrdsc <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(7:12),
    type%in%c(
      "IUUI"#, 
      #"trueX1", "trueX2"
    )
  ) %>% 
  dplyr::select(-s,-unc, -mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k) %>% 
  summarise(
    rr = mean(dsc<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j), 
    k=unique(k)
  ) 
rrdsc


# Rejection rates DM test DSC simulation
rrs_dsc <- 
  dat %>% 
  dplyr::filter(
    setup %in% c(8,10,12),
    type%in%c("pval")
  ) %>% 
  dplyr::select(-unc, -dsc,-mcb, -seed, -i_MC) %>% 
  # be careful !! on which level I have to group?
  group_by(setup,tt,type, k) %>% 
  summarise(
    rr = mean(s<0.1),
    dscinfo = unique(dscinfo), 
    mcbinfo = unique(mcbinfo),
    j=unique(j), 
    k=unique(k)
  ) %>% 
  mutate(
    type="DM",
    panel="dsc"
  )


datall <- 
  bind_rows(
    rrmcb %>% mutate(panel = "mcb", type = "proposed test on MCB difference"),
    rrs_mcb,
    rrdsc %>% mutate(panel = "dsc", type = "proposed test on DSC difference"),
    rrs_dsc
  ) 

lsize <- 1.05

p1<-
  datall %>% 
  # invisible layer to get DSC into legend
  dplyr::filter(panel == "mcb") %>% 
  bind_rows(
    datall %>% 
      dplyr::filter(panel == "mcb") %>% 
      distinct(dscinfo,mcbinfo) %>% 
      mutate(k = 0, rr = Inf, type = "proposed test on DSC difference")
  ) %>% 
  mutate(
    type=factor(type),
    type =recode_factor(
      type,
      "proposed test on MCB difference" = "proposed test on MCB difference",
      "proposed test on DSC difference" = "proposed test on DSC difference",
      "DM test"                     = "DM test",
      .ordered = TRUE   # this fixes the legend order to input order
    )
  ) %>% 
  ggplot(aes(y = rr, x = k, color = type,group = type, linetype=type)) +
  geom_abline(slope = 0, intercept=.1, color ="gray")+
  geom_line(size = lsize) +
  facet_grid(dscinfo ~ mcbinfo, labeller = label_parsed) 
p1

p2 <- 
  datall %>% 
  dplyr::filter(panel == "dsc") %>% 
  # invisible layer to get MCB into legend
  bind_rows(
    datall %>% 
      dplyr::filter(panel == "dsc") %>% 
      distinct(dscinfo,mcbinfo) %>% 
      mutate(k = 0, rr = Inf, type = "proposed test on MCB difference")
  )  %>% 
  mutate(
    type=factor(type),
    type =recode_factor(
      type,
      "proposed test on MCB difference" = "proposed test on MCB difference        ",
      "proposed test on DSC difference" = "proposed test on DSC difference        ",
      "DM test"                     = "DM test        ",
      .ordered = TRUE   # this fixes the legend order to input order
    )
  ) %>% 
  ggplot(aes(y = rr, x = k, color = type,linetype=type)) +
  geom_abline(slope = 0, intercept=.1, color ="gray")+
  geom_line(size = lsize) +
  facet_grid(mcbinfo ~dscinfo , labeller = label_parsed) 
p2


# plot macros
colors <- c("blue", "coral", alpha("black",0.5))
lines <- c("solid", "solid", 12)
labs <- c(
  "proposed test on MCB difference        ",
  "proposed test on DSC difference        ",
  "DM test        "
)

# combine plot
(p1/p2)+
  patchwork::plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme_bw()&
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),     # Remove margins around the legend
    legend.box.margin = margin(-10, 0, 0, 0), 
    legend.key.size = unit(20, "mm"),
    legend.key.width = unit(20, "mm"),
    legend.text = element_text(
      colour = "black",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust=1,
      colour = "black",
      size = 12),
    axis.text.y = element_text(
      colour = "black",
      size = 11
    ),
    strip.text = element_text(
      size = 11
    ),
    axis.title=element_text(
      size=12
    ),
    #  strip.text.x = element_text(margin = margin(.5,.5,.5,.5, "mm")), 
    #   strip.text.y = element_text(margin = margin(.5,.5,.5,.5, "mm"))
  ) &
  scale_y_continuous(breaks=seq(0, 1, by=0.2))&
  # scale_x_continuous(
  #   breaks=seq(0, .14, by=0.02),
  #   limits = c(0, 0.14)
  # )&
  scale_linetype_manual(
    values=lines, 
    labels = labs
  )&
  scale_color_manual(    
    values = colors,
    labels = labs  )&
  xlab(expression(italic(k)))&
  ylab("emprical rejection rate")

ggsave("plots/sim_m.pdf", height=9, width=10)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/sim_m.pdf", 
  height=9,
  width=10
  )


