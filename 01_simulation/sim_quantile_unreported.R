rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")

# k_=0
# alpha=0.01
# s_=1
# tt=1000


#K <- seq(0.0,.5, length.out=11)
K <- seq(0.0,3,.1)
core.max <- 8
cl <- makeCluster(min(parallel::detectCores() - 1, core.max))
registerDoParallel(cl)
TT <- 2^(9:11)
ALPHA <- c(0.01,0.1)
MCreps <- 2500
t0 <- Sys.time()
MCsim <- foreach(
  i_MC = 1:MCreps,
  #.errorhandling = "pass",
  .packages = c("dplyr", "tibble", "tidyr", "SDI", "purrr")
) %dopar% {
  
  set.seed(i_MC)
  # load functions
  source("../00_functions/funs_simulation.R")
  
  # loop over sample size
  out_tt <- NULL
  for (alpha in ALPHA) {
    for(tt in TT){
      
      # set SE loss and respective identification function
      scoreingfunction = function(x,y, alpha) (1 * (y < x) - alpha) * (x - y)
      identificationfunction = function(x,y, alpha) (1 * (y < x) - alpha)
      
      # loop over setups and k
      out_ks <- NULL
      for(k_ in K){
        for (s_ in c(3)) {
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
          
          
          # obtain dec values an variances
          sdi<- 
            qSDI(
              X_1 = dat$qX1, 
              X_2 = dat$qX2, 
              Y = dat$Y, 
              S=scoreingfunction,
              V=identificationfunction, 
              alpha = alpha, 
              vcov_estimator = sandwich::vcovHAC, 
              tol = 40
            )
          
          
          dec1 <- sdi$asyvardm$dec_1
          dec2 <- sdi$asyvardm$dec_2
          
          
          vars <-  sdi$asyvardm$asy_vars
          # obtain p-values
          pvals <- sdi$asyvardm$pvals
          
          j=NA # add from simulation in large data
          
          
          out <-
            bind_cols(vars, pvals) %>%
            pivot_longer(
              cols = everything(),
              names_to = c(".value", "type"),
              names_sep = "_"
            ) %>%
            bind_rows(
              true, 
              tibble(
                type = "IUUI", 
                mcb = sdi$mcbIUUI,  
                dsc = sdi$dscIUUI
              ), 
              tibble(
                type = "pval_QF_X1", # quadratische Form forecast 1
                s = NA, mcb=sdi$mcbnulltest_X1$Wald.pval.MCB, dsc=sdi$dscnulltest_X1$Wald.pval.DSC
              ), 
              tibble(
                type = "pval_QF_X2", # quadratische Form forecast 2
                s = NA, mcb=sdi$mcbnulltest_X2$Wald.pval.MCB, dsc=sdi$dscnulltest_X2$Wald.pval.DSC
              ), 
              tibble(
                type= "bonf_pval", 
                s = NA, 
                mcb= 2*pmin(sdi$mcbnulltest_X1$Wald.pval.MCB, sdi$mcbnulltest_X2$Wald.pval.MCB),
                dsc= 2*pmin(sdi$dscnulltest_X1$Wald.pval.DSC, sdi$dscnulltest_X2$Wald.pval.DSC)
              )
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
              i_MC = i_MC,
              seed = i_MC, #+k_+s_
              j= j,
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
  }# close loop over alpha
  
  out_tt
}
t1 <- Sys.time()
t1-t0# 4.4h,  5000 reps 
stopCluster(cl)
dat <- do.call(rbind, MCsim)
saveRDS(dat,"data/sim_q_parameterized_size_details.rds") 

# Plot size detail  -----
dat <- readRDS("data/sim_q_parameterized_size_details.rds")
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
  # be careful !! on which level to group?
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
    #rrmcb  %>% mutate(panel = "mcb", type = "proposed test on MCB difference"),
    #rrs_mcb %>% mutate(panel = "mcb", type = "DM"),
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
  mutate(
    alpha_lab = paste0("alpha == ", alpha), 
    tt_lab = paste0("T == ", tt)
  ) 
colors <- c("blue", "coral", alpha("black",0.5))
lines <- c("solid", "solid", "dotted")
datall %>% 
  ggplot(aes(y = rr, x = k, color = type,linetype = type)) +
  geom_abline(slope = 0, intercept = .1, color = "gray") +
  geom_line() +
  facet_grid(alpha_lab ~ tt + mcbinfo + dscinfo, labeller = label_parsed) +
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
  ) +
  xlab(expression(italic(k))) +
  ylab("empirical rejection rate")+
  scale_linetype_manual(
    values=lines[-1]
  )+
  scale_color_manual(    
    values = colors[-1]
  )
