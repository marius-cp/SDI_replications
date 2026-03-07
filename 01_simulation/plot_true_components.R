rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")

# Globals ----

ord <- c("MCB1","DSC1","MCB2","DSC2")

labs_math <- c(
  MCB1 = expression(bar(MCB)[1]),
  DSC1 = expression(bar(DSC)[1]),
  MCB2 = expression(bar(MCB)[2]),
  DSC2 = expression(bar(DSC)[2])
)

cols_forecasters <- c(
  MCB1 = "blue",
  DSC1 = "magenta",
  MCB2 = "coral",
  DSC2 = "black"
)


lt_map <- c(MCB1 = "solid", DSC1 = "solid", MCB2 = "dotted", DSC2 = "dotted")

# True Mean ----
K <- seq(0.0,.5, length.out=11)
core.max <- 1
cl <- makeCluster(min(parallel::detectCores() - 1, core.max))
registerDoParallel(cl)
TT <- 500 
MCreps <- 1
t0 <- Sys.time()
MCsim <- foreach(
  i_MC = 1:MCreps,
  #.errorhandling = "pass",
  .packages = c("dplyr", "tibble", "tidyr", "SDI")
) %dopar% {
  
  set.seed(i_MC)
  # load functions
  source("../00_functions/funs_simulation.R")
  
  # loop over sample size
  out_tt <- NULL
  for(tt in TT){
    
    # set SE loss and respective identification function
    scoreingfunction = function(x,y) (x-y)^2
    identificationfunction = function(x,y) 2*(x-y)
    spp = function(x,y) 2
    
    # loop over setups and k
    out_ks <- NULL
    for(k_ in K){
      for (s_ in 1:12) {
        #set.seed(i_MC)
        
        param <- parameters_simset_m(s_=s_,k_=k_)
        
        dat <-
          dgp(
            t = tt,
            gamma = param$gamma,
            beta = param$beta,
            delta = param$delta,
            xi = param$xi
          ) 
        
        trueX1 <- se_decomp_prop_sim(
          gamma = param$gamma,
          beta = param$beta, 
          parameterforecaster = param$delta
        )
        trueX2 <- se_decomp_prop_sim(
          gamma = param$gamma,
          beta = param$beta, 
          parameterforecaster = param$xi
        )
        true_components <- bind_rows(
          data.frame(trueX1) %>% mutate(type="true_X1"),
          data.frame(trueX2) %>% mutate(type="true_X2")
        ) %>% 
          select(-varsigma) %>% 
          mutate(s=mcb-dsc+unc)
        
     
        ifelse(
          s_ <= 6, 
          # MCB simulation j 
          j <-  true_components$mcb[2],
          # DSC simulation j 
          j <-  true_components$dsc[2]
        )
        
        out <-
          true_components %>%
          pivot_longer(
            cols = c("mcb", "dsc", "unc", "s"),
          ) %>% 
          mutate(
            tt = tt,
            setup = s_, # remember that s is already defined as score!!!
            k=k_,
            i_MC = i_MC,
            seed = i_MC, #+k_+s_
            j= j,
            mcbinfo = factor(1, levels=1, labels=param$expression_mcb), 
            dscinfo = factor(1, levels=1, labels=param$expression_dsc)
          ) 
        out_ks <- rbind(out_ks,out)
        
      } # close loop over k
    }# close loop s
    out_tt <- rbind(out_tt,out_ks)
  }# close loop over sample size T
  out_tt
}
t1 <- Sys.time()
t1-t0
stopCluster(cl)


dat <- do.call(rbind, MCsim)

out_long_true_mcb <-
  dat %>%
  filter(type %in% c("true_X1","true_X2")) %>% 
  filter(name != "unc", name != "s", setup <= 6) %>%
  mutate(
    series = case_when(
      name == "mcb" & type == "true_X1" ~ "MCB1",
      name == "mcb" & type == "true_X2" ~ "MCB2",
      name == "dsc" & type == "true_X1" ~ "DSC1",
      name == "dsc" & type == "true_X2" ~ "DSC2",
      TRUE ~ NA_character_
    ),
    series = factor(series, levels = c("MCB1", "DSC1", "MCB2", "DSC2"))
  ) %>%
  arrange(mcbinfo, dscinfo, k, series)

pmmcb <- 
ggplot(out_long_true_mcb, aes(x = k, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  #facet_nested(alpha ~ mcbinfo + dscinfo, labeller = label_parsed, scales = "free") +
  facet_grid2(
    dscinfo ~ mcbinfo,
    labeller = label_parsed,
   # scales = "free_y",
    #independent = "y"
  )+
  scale_color_manual(
    breaks = ord,
    values = cols_forecasters[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  scale_linetype_manual(
    breaks = ord,
    values = lt_map[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  guides(
    linetype = "none",
    color = guide_legend(
      title = NULL,
      override.aes = list(linetype = unname(lt_map[ord]))
    )
  ) +
  theme_bw() +
  xlab(expression(italic(k))) 

out_long_true_dsc <-
  dat %>%
  filter(type %in% c("true_X1","true_X2")) %>% 
  filter(name != "unc", name != "s", setup > 6) %>%
  mutate(
    series = case_when(
      name == "mcb" & type == "true_X1" ~ "MCB1",
      name == "mcb" & type == "true_X2" ~ "MCB2",
      name == "dsc" & type == "true_X1" ~ "DSC1",
      name == "dsc" & type == "true_X2" ~ "DSC2",
      TRUE ~ NA_character_
    ),
    series = factor(series, levels = c("MCB1", "DSC1", "MCB2", "DSC2"))
  ) %>%
  arrange(mcbinfo, dscinfo, k, series)

pmdsc <- 
ggplot(out_long_true_dsc, aes(x = k, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  #facet_nested(alpha ~ mcbinfo + dscinfo, labeller = label_parsed, scales = "free") +
  facet_grid2(
    mcbinfo ~ dscinfo,
    labeller = label_parsed,
    #scales = "free_y",
    #independent = "y"
  )+
  scale_color_manual(
    breaks = ord,
    values = cols_forecasters[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  scale_linetype_manual(
    breaks = ord,
    values = lt_map[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  guides(
    linetype = "none",
    color = guide_legend(
      title = NULL,
      override.aes = list(linetype = unname(lt_map[ord]))
    )
  ) +
  theme_bw() +
  xlab(expression(italic(k))) 


(pmmcb/pmdsc)+
  patchwork::plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(12, "mm"),
    legend.key.size = unit(12, "mm"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )&
  ylab("Population value")

ggsave("plots/m_true.pdf", height=9, width=10)
ggsave("/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/m_true.pdf", height=9, width=10)




# True Quantile ----
K <- seq(0.0,1, length.out=21)
core.max <- 5
cl <- makeCluster(min(parallel::detectCores() - 1, core.max))
registerDoParallel(cl)
TT <- 100000
ALPHA <- c(0.01,0.05,0.1,0.25,.5)
MCreps <- 8
t0 <- Sys.time()
compare <- foreach(
  al = 1:length(ALPHA),
  #.errorhandling = "pass",
  .packages = c("dplyr", "tibble", "tidyr", "SDI", "purrr")
) %dopar% {
  seed<-123
  set.seed(seed = seed)
  out_tt <- NULL
  alpha <- ALPHA[al]
  for(tt in TT){ # we can use large samples here to double check if true MCB calc matchs the emprical value
    # set SE loss and respective identification function
    scoreingfunction = function(x,y, alpha) (1 * (y < x) - alpha) * (x - y)
    identificationfunction = function(x,y, alpha) (1 * (y < x) - alpha)
    
    # loop over setups and k
    out_ks <- NULL
    for(k_ in K){
      for (s_ in seq(1,4)) {
        #set.seed(i_MC)
        param <- parameters_simset_q(s_=s_, k_=k_)
        
        # update xi0 if we are in DGP 4
        ifelse(
          s_==4,
          param$xi[1] <- getxi0(ALPHA=alpha, K=k_),
          param$xi[1] <- param$xi[1] 
        )
        
        dat <-
          dgp(
            t = tt,
            gamma = param$gamma,
            beta = param$beta,
            delta = param$delta,
            xi = param$xi, 
            alpha = alpha
          ) 
        
        true1 <- mcb_dsc_checkloss_closed(
          alpha  = alpha,
          beta   = unique(param$beta),
          gamma  = param$gamma,
          delta  = param$delta[-1],
          delta0 = param$delta[1]
        )
        true2 <- mcb_dsc_checkloss_closed(
          alpha  = alpha,
          beta   = unique(param$beta),
          gamma  = param$gamma,
          delta  = param$xi[-1],
          delta0 = param$xi[1]
        )
        
        true <- bind_rows(
          as_tibble(true1[c("mcb", "dsc", "unc")]) %>% mutate(type = "true1"),
          as_tibble(true2[c("mcb", "dsc", "unc")]) %>% mutate(type = "true2")
        )
        
        estimate1 <-
          qdecomposition(
            S = scoreingfunction,
            alpha=alpha,
            x=dat$qX1,
            y=dat$Y, 
            x_rc = as.numeric(fitted.values(quantreg::rq(Y ~ qX1, tau = alpha, data = dat)))
          )
        
        estimate2 <-
          qdecomposition(
            S = scoreingfunction,
            alpha=alpha,
            x=dat$qX2,
            y=dat$Y, 
            x_rc = as.numeric(fitted.values(quantreg::rq(Y ~ qX2, tau = alpha, data = dat)))
          )
        
        estimate <- bind_rows(
          as_tibble(estimate1[c("mcb", "dsc", "unc")]) %>% mutate(type = "estimate1"),
          as_tibble(estimate2[c("mcb", "dsc", "unc")]) %>% mutate(type = "estimate2")
        )
        estimate
        
        differences1 <- 
          true %>%
          filter(type == "true1") %>%
          select(mcb, dsc, unc) %>%
          rename_with(~ paste0(.x, "_true")) %>%
          bind_cols(
            estimate %>%
              filter(type == "estimate1") %>%
              select(mcb, dsc, unc) %>%
              rename_with(~ paste0(.x, "_est"))
          ) %>%
          # mutate(
          #   mcb = mcb_est - mcb_true,
          #   dsc = dsc_est - dsc_true,
          #   unc = unc_est - unc_true, 
          #   type = "diff"
          # ) %>% 
          mutate(
            mcb = 100 * (mcb_est - mcb_true) / mcb_true,
            dsc = 100 * (dsc_est - dsc_true) / dsc_true,
            unc = 100 * (unc_est - unc_true) / unc_true,
            type = "diff1"
          ) %>% 
          select(mcb:type)
        
        differences2 <- 
          true %>%
          filter(type == "true2") %>%
          select(mcb, dsc, unc) %>%
          rename_with(~ paste0(.x, "_true")) %>%
          bind_cols(
            estimate %>%
              filter(type == "estimate2") %>%
              select(mcb, dsc, unc) %>%
              rename_with(~ paste0(.x, "_est"))
          ) %>%
          # mutate(
          #   mcb = mcb_est - mcb_true,
          #   dsc = dsc_est - dsc_true,
          #   unc = unc_est - unc_true, 
          #   type = "diff"
          # ) %>% 
          mutate(
            mcb = 100 * (mcb_est - mcb_true) / mcb_true,
            dsc = 100 * (dsc_est - dsc_true) / dsc_true,
            unc = 100 * (unc_est - unc_true) / unc_true,
            type = "diff2"
          ) %>% 
          select(mcb:type)
        
        differences <- bind_rows(differences1,differences2)
        
        out <-
          bind_rows(
            true, 
            estimate,
            differences
          ) %>% 
          mutate(
            type = case_when(
              type == "1" ~ "X1",
              type == "2" ~ "X2",
              TRUE ~ type
            ), 
            tt = tt,
            setup = s_, # remember that s is already defined as score!!!
            k=k_,
            seed = seed, #+k_+s_
            alpha = alpha,
            mcbinfo = factor(1, levels=1, labels=param$expression_mcb), 
            dscinfo = factor(1, levels=1, labels=param$expression_dsc), 
            xi0=param$xi[1]
          ) 
        out_ks <- rbind(out_ks,out)
        
      } # close loop over k
    }# close loop s
    out_tt <- rbind(out_tt,out_ks)
  }# close loop over sample size T
  out_tt
}# close loop over alpha
t1 <- Sys.time()
t1-t0
stopCluster(cl)
dat <- do.call(rbind, compare)
saveRDS(dat,"data/sim_q_compare.rds") 

out_tt <-readRDS("data/sim_q_compare.rds")

## plot ----
out_long_true <-
  out_tt %>%
  filter(type %in% c("true1","true2")) %>% 
  select(-unc) %>%
  pivot_longer(cols = mcb:dsc, names_to = "name", values_to = "value") %>%
  mutate(
    series = case_when(
      name == "mcb" & type == "true1" ~ "MCB1",
      name == "mcb" & type == "true2" ~ "MCB2",
      name == "dsc" & type == "true1" ~ "DSC1",
      name == "dsc" & type == "true2" ~ "DSC2",
      TRUE ~ NA_character_
    ),
    series = factor(series, levels = c("MCB1", "DSC1", "MCB2", "DSC2"))
  ) %>%
  arrange(alpha, mcbinfo, dscinfo, k, series) %>% 
  mutate(alpha_lab = paste0("alpha == ", alpha))

ggplot(out_long_true, aes(x = k, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.8) +
  #facet_nested(alpha ~ mcbinfo + dscinfo, labeller = label_parsed, scales = "free") +
  facet_grid2(
    alpha_lab ~ mcbinfo + dscinfo,
    labeller = label_parsed,
    scales = "free_y",
    independent = "y"
  )+
  scale_color_manual(
    breaks = ord,
    values = cols_forecasters[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  scale_linetype_manual(
    breaks = ord,
    values = lt_map[ord],
    labels = labs_math[ord],
    name   = NULL
  ) +
  guides(
    linetype = "none",
    color = guide_legend(
      title = NULL,
      override.aes = list(linetype = unname(lt_map[ord]))
    )
  ) +
  theme_bw() +
  xlab(expression(italic(k))) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(12, "mm"),
    legend.key.size = unit(12, "mm"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )+
  ylab("Population value")
ggsave("plots/q_true.pdf", height = 11, width = 11)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/q_true.pdf", 
  height=11, 
  width=11
  )


out_tt %>% 
  filter(type %in% c("diff1", "diff2")) %>% 
  select(-unc) %>%
  pivot_longer(cols = mcb:dsc, names_to = "name", values_to = "value") %>% 
  mutate(
    alpha_lab = paste0("alpha == ", alpha),
    type=ifelse(type=="diff1", "X1", "X2")
  ) %>%
  mutate(
    series = case_when(
      name == "mcb"  ~ "MCB",
      name == "dsc" ~ "DSC",
      TRUE ~ NA_character_
    ),
    series = factor(series, levels = c("MCB", "DSC"))
  ) %>% 
  ggplot(
    aes(x = k, y = value, color = series)
  )+
  geom_abline(intercept = 0, slope = 0)+
  facet_nested(
    type + alpha_lab ~ mcbinfo + dscinfo,
    labeller = label_parsed,
    scales = "free_y",
    independent = "y"
  )+
  geom_line(linewidth = 0.8) +
  theme_bw() +
  xlab(expression(italic(k))) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(12, "mm"),
    legend.key.size = unit(12, "mm"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )+
  ylab("percentage difference (estimate − true) / true [%]")

## plot constants ----
readRDS("data/sim_q_compare.rds") %>% 
  filter(type=="true2", setup==4) %>% 
  mutate(alpha_lab = paste0("alpha == ", alpha)) %>% 
  ggplot(aes(x = k, y = xi0)) +
  geom_line(linewidth = 0.8) +
  facet_nested(
    #    .~ mcbinfo + dscinfo+ alpha_lab ,
    .~alpha_lab ,
    labeller = label_parsed,
    scales = "free_y",
    independent = "y"
  )+
  ylab(expression(xi[0]))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y   = unit(1, "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(12, "mm"),
    legend.key.size = unit(12, "mm"),
    legend.text = element_text(colour = "black", size = 12),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 11),
    strip.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )+
  xlab(expression(italic(k))) 
ggsave("plots/q_xi0_choice.pdf", height = 2.5, width = 10)
ggsave(
  "/Users/mp/Library/CloudStorage/Dropbox/Apps/Overleaf/Statistical Inference for Score Decompositions/fig/q_xi0_choice.pdf", 
  height = 2.5, width = 10
)

