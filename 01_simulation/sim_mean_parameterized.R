rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(SDI)
library(patchwork)
library(doParallel)
source("../00_functions/funs_simulation.R")

# MAIN Simulation ---- 
K <- seq(0.0,.5, length.out=11)
core.max <- 8
cl <- makeCluster(min(parallel::detectCores() - 1, core.max))
registerDoParallel(cl)
TT <- 500 
MCreps <- 5000
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
        
        # obtain dec values an variances
        sdi<- 
          SDI(
            X_1 = dat$X1, 
            X_2 = dat$X2, 
            Y = dat$Y, 
            S=scoreingfunction,
            V=identificationfunction, 
            Spp = spp, 
            vcov_estimator = sandwich::vcovHAC, 
            tol = 40
          )
        
        mz1 <- lm(Y~X1, data = dat)
        mz2 <- lm(Y~X2, data = dat)
        
        dec1 <- sdi$asyvardm$dec_1
        dec2 <- sdi$asyvardm$dec_2
        
        
        vars <-  sdi$asyvardm$asy_vars
        # obtain p-values
        pvals <- sdi$asyvardm$pvals
        
        
        ifelse(
          s_ <= 6, 
          # MCB simulation j 
          j <-  true_components$mcb[2],
          # DSC simulation j 
          j <-  true_components$dsc[2]
        )
        
        
        bind_cols(vars, pvals) %>%
          pivot_longer(
            cols = everything(),
            names_to = c(".value", "type"),
            names_sep = "_"
          ) %>%
          bind_rows(
            true_components)
        
        out <-
          bind_cols(vars, pvals) %>%
          pivot_longer(
            cols = everything(),
            names_to = c(".value", "type"),
            names_sep = "_"
          ) %>%
          bind_rows(
            true_components, 
            tibble(
              type = "IUUI", 
              mcb = sdi$mcbIUUI,  
              dsc = sdi$dscIUUI
            ), 
            tibble(
              type = "pval_QF_X1", # quadratische Form forecast 1
              s = NA, mcb=sdi$mcbnulltest_X1$mcb_pval, dsc=sdi$dscnulltest_X1$dsc_pval
            ), 
            tibble(
              type = "pval_QF_X2", # quadratische Form forecast 2
              s = NA, mcb=sdi$mcbnulltest_X2$mcb_pval, dsc=sdi$dscnulltest_X2$dsc_pval
            ), 
            tibble(
              type= "bonf_pval", 
              s = NA, 
              mcb= 2*pmin(sdi$mcbnulltest_X1$mcb_pval, sdi$mcbnulltest_X2$mcb_pval),
              dsc= 2*pmin(sdi$dscnulltest_X1$dsc_pval, sdi$dscnulltest_X2$dsc_pval)
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
t1-t0# 6h, 5000 reps
stopCluster(cl)


dat <- do.call(rbind, MCsim)
saveRDS(dat,"data/sim_m_parameterized.rds") 

