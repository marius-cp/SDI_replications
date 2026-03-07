rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")

# sim----
#K <- seq(0.0,.5, length.out=11)
K <- seq(0.0,1, length.out=21)
core.max <- 8
cl <- makeCluster(min(parallel::detectCores() - 1, core.max))
registerDoParallel(cl)
TT <- 500
ALPHA <- c(0.01,0.05,0.1,0.25,.5)
MCreps <- 5000
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
saveRDS(dat,"data/sim_q_parameterized.rds") 

# Split large simulation output file for GitHub compatibility
#
# The simulation output `sim_q_parameterized.rds` contains ~21 million rows and
# results in a file size (~300 MB) that exceeds GitHub's recommended file size
# limits. To allow the data to be included in the repository, the file is split
# into four smaller parts.
dat <- readRDS("data/sim_q_parameterized.rds")
n <- nrow(dat)
idx <- split(seq_len(n), cut(seq_len(n), 4, labels = FALSE))
for(i in seq_along(idx)) {
  saveRDS(dat[idx[[i]], ], paste0("data/sim_q_parameterized_part", i, ".rds"))
}
