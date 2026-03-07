# Backtest functions ----
#
# This script contains custom backtesting functions for forecast evaluation.
# Parts of the implementation are inspired by or adapted from existing open-source
# code, in particular:
#
# - rugarch: https://rdrr.io/cran/rugarch/src/R/rugarch-tests.R
# - Nolde et al. (Elicitability and Backtesting):
#   https://github.com/nnolde/Elicitability-and-Backtesting
#


trafficlight <- function(y,VaR,VaR_level){
  y=tail(y,length(VaR))
  n=length(y)
  
  hit=numeric(n)
  hit[y < VaR] = 1
  no_fail=sum(hit)
  
  prob = pbinom(q=no_fail,size = n,prob = 1-VaR_level)
  
  typeIprob=1-prob+dbinom(x=no_fail,size = n,prob = 1-VaR_level)
  
  if (prob <= 0.95) {
    TL_color="green"
  } else if (prob>0.95 && prob <= 0.9999){
    TL_color="yellow"
  } else TL_color="red"
  
  res = list()
  res$prob = prob
  res$typeI = typeIprob
  res$color = TL_color
  res$absHits = no_fail
  res$relHits = no_fail/n
  return(res)
}

# now with HAC for QR!!
VQRbacktest_new <- function(alpha = 0.05, actual, VaR,
                            qr_cov = c("default","hac"),
                            ker = 3L,
                            bandwidth = NULL){
  
   VQR <- SDI::MZ.test.quantile(
     x=VaR, y=actual, alpha = alpha, 
     cov = qr_cov, 
     ker = 3L, 
     bandwidth = NULL
     )
  
  list(
    alpha = alpha,
    qr_cov = qr_cov,
    coef = VQR$estimate,
    stat = as.numeric(VQR$WaldStat.MCB),
    pval = as.numeric(VQR$Wald.pval.MCB)
    )
}



VQRbacktest = function(alpha = 0.05, actual, VaR, test = "original"){
  
  ans = list()
  
  if(test=="original"){
    tmp = SDI::MZ.test.quantile(alpha = alpha, y = actual, x  = VaR)
    ans$pval = tmp$Wald.pval.MCB
    ans$fitvals = tmp$fittedvalues
    ans$estimate = tmp$estimate
  }
  
  if(test == "wald"){
    r = VaR
    x=actual
    
    lev = 1-alpha
    
    
    n = length(x) # out-of-sample size
    # setting up the identification function
    v = (r<=x) - lev 
    
    N=n 
    vt=v #v[t], t=2,...,n, identification function used in test statistic computation
    h=rbind(rep(1,N),r)
    
    q=dim(h)[1]
    
    zbar=(h %*% (vt))/N # q x 1 , here vt is shifted one step into future
    Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
    for (t in 1:dim(h)[2])
    {			
      Omega = Omega + vt[t]^2*(h[,t]%*%t(h[,t]))	
      
    }
    Omega = Omega/N
    Omega.inv = solve(Omega)
    
    tn.cond <- N*t(zbar)%*% Omega.inv %*% zbar 	#test statistic
    pv.cond <- pchisq(tn.cond, df=q, lower.tail=FALSE) # p-value	
    
    ans$crit = tn.cond
    ans$pval = pv.cond
    
  }
  
  
  return(ans)
}



regression_based_backtest_ols <- function(
    data,
    formulas = list(
      UC  = hit_demeaned ~ 1,
      CC  = hit_demeaned ~ 1 + hit_lag1,
      DQ  = hit_demeaned ~ 1 + value + hit_lag1 + hit_lag2 + hit_lag3 + hit_lag4,
      DQv = hit_demeaned ~ 1 + hit_lag1 + hit_lag2 + hit_lag3 + hit_lag4
    ),
    vcov_fun = sandwich::vcovHAC,
    test = c("F", "Chisq"),
    return_full = FALSE
) {
  test <- match.arg(test)
  
  run_one <- function(formula) {
    reg <- stats::lm(formula, data = data)
    V   <- vcov_fun(reg)
    
    cn    <- names(stats::coef(reg))
    restr <- paste0(cn, " = 0")
    
    lh <- car::linearHypothesis(reg, restr, vcov. = V, test = test)
    
    pval <- lh[2, "Pr(>F)"]
    if (is.na(pval) && "Pr(>Chisq)" %in% colnames(lh)) {
      pval <- lh[2, "Pr(>Chisq)"]
    }
    
    if (return_full) {
      list(
        p_value = as.numeric(pval),
        test_table = lh,
        reg = reg,
        vcov = V
      )
    } else {
      as.numeric(pval)
    }
  }
  
  res <- lapply(formulas, run_one)
  
  if (return_full) {
    return(res)
  }
  
  # one-row data.frame with named columns
  as.data.frame(res, row.names = "p_value")
}

backtesttab_wrapper <- function(data, alpha) {
  
  stopifnot(all(c("ret", "value") %in% names(data)))
  
  ## --- regression-based backtests (OLS / HAC) ---
  olsregs <- regression_based_backtest_ols(data)
  
  ## --- VQR backtest ---
  # VQR <- VQRbacktest(
  #   alpha  = unique(data$alpha),
  #   actual = data$ret,
  #   VaR    = data$value,
  #   test   = "original"
  # )
  VQR <- VQRbacktest_new(unique(data$alpha), data$ret, data$value, qr_cov=c("hac"))
  
  
  ## --- Basel traffic light ---
  TL <- trafficlight(
    y         = data$ret,
    VaR       = data$value,
    VaR_level = 1 - unique(data$alpha)
  )
  
  NZ <- NZtest(x=-data$ret,r=-data$value, lev=1-unique(data$alpha), hac=TRUE)
  
  ## --- assemble output ---
  tibble::tibble(
    UC  = olsregs$UC,
    CC  = olsregs$CC,
    DQ  = olsregs$DQ,
    DQv = olsregs$DQv,
    VQR    = as.numeric(VQR$pval),
    Basel  = 1 - TL$prob,
    Baselc = TL$color,
    NZ = NZ$pv.avg,
    relHits= TL$relHits*100 # in percent
  )
}

# ===================================================================
# ONE-SIDED TEST of super-calibration "cct.1s.VaR"
# test statistic tau2 for average tests and tau4 for conditional tests
# ===================================================================

NZtest <- function(
    x, # rel
    r, # forecast
    lev=0.95, hac=FALSE) # super-calbration
{
  n = length(x) # out-of-sample size
  
  # setting up appropriate identification functions for the risk measure
  v = (r>=x) - lev
  

  NZ_HACestimator <- function(x)
  {
    n = length(x)
    m=ceiling(2*sqrt(n))
    gam = acf(x,lag.max=m,type="covariance",plot=F)$acf
    k1 = 1:ceiling(m/2)
    k2 = (ceiling(m/2)+1):m
    
    lam = c(1, 2*(1-6*(k1/m)^2+6*(k1/m)^3),2*2*(1-k2/m)^3)
    (gam %*% lam)
  }
  
  # average calibration test
  if(hac) Omega = NZ_HACestimator(v)
  else{ Omega = (v %*% v)/n  # zero truncation lag
  }
  
  tn.avg <- sqrt(n)*mean(v)/sqrt(Omega) 	#test statistic
  pv.avg <- pnorm(tn.avg) # p-value	
  
  # conditional calibration tests	
  
  ## same test function as for the two-sided test
  N=n
  vt=v
  #h=rbind(rep(1,N),abs(r))
  h=rbind(rep(1,N),(r))
  
  q=dim(h)[1]
  
  zbar=(h %*% (vt))/N # q x 1
  Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
  for (t in 1:N)
  {		
    Omega = Omega + vt[t]^2*(h[,t]%*%t(h[,t]))	# V is a scalar for k=1
  }
  Omega = Omega/N
  
  # Hommel's (1983) procedure
  tn <- sqrt(N) * diag(Omega)^(-1/2) * zbar 	#test statistic
  pi <- sort(pnorm(tn))
  cq <- sum(1/(1:q))
  pv.cond <- min(q*cq*min(pi/(1:q)),1)
  
  pv.condB <- min(q*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
  
  return(list(pv.avg = pv.avg, pv.cond = pv.cond, pv.condB = pv.condB))
}

