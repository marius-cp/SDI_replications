#' @title Data Generating Process (DGP)
#'
#' @description
#' Simulates a simple quantile data generating process for \eqn{Y_t} with two
#' forecasters \eqn{X_{1t}} and \eqn{X_{2t}} based on AR(1) predictors
#' \eqn{(K_t, L_t, M_t,  N_t)}. Returns both quantile forecasts based on the  mean forecasts using the location-scale property. 
#'
#' @param t Sample size. Default is 1000.
#' @param gamma Coefficients of predictors in the DGP of \eqn{Y_t}.
#' @param beta AR(1) parameters for \eqn{(M,K,L,N)}.
#' @param delta Coefficients for forecaster 1.
#' @param xi Coefficients for forecaster 2.
#' @param alpha Quantile level (e.g., 0.25).
#'
#' @return
#' A \code{data.frame} with columns \code{Y}, \code{X1}, \code{X2},
#' \code{qX1}, and \code{qX2}.
#'
#' @examples
#' set.seed(123)
#' qdgp(t = 500, alpha = 0.1)
#'
dgp <- function(
    t = 1000,
    gamma = c(1,1,1, 1),
    beta = c(rep(0.25,4)), 
    delta = c(0,0,0,1,0), # delta_0, and rest \bm{\delta}
    xi = c(0,0,1,0,0.9,0), # xi_0, and rest \bm{\xi}
    alpha=0.25
){
  
  # series K
  K <- arima.sim(list(order=c(1,0,0), ar=beta[2]), n=t)
  L <- arima.sim(list(order=c(1,0,0), ar=beta[3]), n=t)
  M <- arima.sim(list(order=c(1,0,0), ar=beta[1]), n=t)
  N <- arima.sim(list(order=c(1,0,0), ar=beta[3]), n=t)
  # stack into V_t = (K_t, L_t, M_t, N_t)'
  V <- cbind(K, L, M, N)
  # innovation for Y (you can change sd if you want)
  u <- rnorm(t)
  # generate Y_t = V_t' gamma + u
  Y <- as.numeric(V %*% gamma + u)
  # fcast 1
  X1 <- delta[1] + V %*%delta[2:5]
  # fcast 2
  X2 <- xi[1] + V %*%xi[2:5]
  # quantile forecast
  zalpha <- qnorm(alpha)        # if F = N(0,1)
  sigma  <- 1 
  qX1 <- X1 + zalpha * sigma
  qX2 <- X2 + zalpha * sigma
  
  out <- data.frame(
    Y, X1, X2, qX1, qX2
  )
  
  return(out)
}

# sets parameters for mean simulation
parameters_simset_m <- function(s_,k_){
  # MCB simulations:
  # A1a
  if(s_==1){
    beta_=c(rep(0.25,4))
    gamma_=c(0,1/4,1/4,0)
    delta_=c(0,0,0,k_/sqrt(2)+1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    #expression_mcb <- expression(bar(MCB)[1] == 0 ~~ plain("and") ~~ bar(MCB)[2] == frac(8 * k^2, 15))
    expression_mcb <- expression({bar(MCB)[1](italic(k)) == bar(MCB)[2](italic(k))} == 8 * italic(k)^2 /15)
    expression_dsc <- expression(0 < {bar(DSC)[1] < bar(DSC)[2]})
  }
  # A1b
  if(s_==2){
    beta_=c(rep(0.25,4))
    gamma_=c(1/4,1/4,1/4,0)
    delta_=c(0,k_/2+1/4,0,k_/2+1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    expression_mcb <- expression({bar(MCB)[1](italic(k)) == bar(MCB)[2](italic(k))} == 8 * italic(k)^2 /15)
    expression_dsc <- expression(0 < {bar(DSC)[1] == bar(DSC)[2]})
  }
  # A2a
  if(s_==3){
    beta_=c(rep(0.25,4))
    gamma_=c(0,1/4,1/4,0)
    delta_=c(0,0,0,1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    expression_mcb = expression({{bar(MCB)[1] == 0} ~~ plain(and) ~~bar(MCB)[2](italic(k)) == 8 * italic(k)^2 /15})
    expression_dsc <- expression(0 < {bar(DSC)[1] < bar(DSC)[2]})
    
  }
  # A2b
  if(s_==4){
    beta_=c(rep(0.25,4))
    gamma_=c(1/4,1/4,1/4,0)
    delta_=c(0,1/4,0,1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    expression_mcb = expression({{bar(MCB)[1] == 0} ~~ plain(and) ~~bar(MCB)[2](italic(k)) == 8 * italic(k)^2 /15})
    expression_dsc <- expression(0 < {bar(DSC)[1] == bar(DSC)[2]})
  }
  # A3a
  if(s_==5){
    beta_=c(rep(0.25,4))
    gamma_=c(0,1/4,1/4,0)
    delta_=c(0,0,0,k_+1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    expression_dsc <- expression(0 < {bar(DSC)[1] < bar(DSC)[2]})
    expression_mcb <- expression({bar(MCB)[1](italic(k)) == 2*bar(MCB)[2](italic(k))} == 16 * italic(k)^2 /15)
    
  }
  # A3b
  if(s_==6){
    beta_=c(rep(0.25,4))
    gamma_=c(1/4,1/4,1/4,0)
    delta_=c(0,k_/sqrt(2)+1/4,0,k_/sqrt(2)+1/4,0)
    xi_=c(0,0,k_/2+1/4,k_/2+1/4,0)
    expression_dsc <- expression(0 < {bar(DSC)[1] == bar(DSC)[2]})
    expression_mcb <- expression({bar(MCB)[1](italic(k)) == 2*bar(MCB)[2](italic(k))} == 16 * italic(k)^2 /15)
  }
  
  # DSC simulations:
  # B1a
  if(s_==7){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(0,k_/2+1/(4*sqrt(2)),0,0,k_/2+1/(4*sqrt(2)))
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    expression_dsc <- expression({bar(DSC)[1](italic(k)) == bar(DSC)[2](italic(k))} == 8 * italic(k)^2 /15)
    expression_mcb <- expression(0 < {bar(MCB)[1] < bar(MCB)[2]})
  }
  # B1b
  if(s_==8){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(0,k_/2+1/4,0,0,k_/2+1/4)
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    expression_dsc <- expression({bar(DSC)[1](italic(k)) == bar(DSC)[2](italic(k))} == 8 * italic(k)^2 /15)
    expression_mcb <- expression(0 < {bar(MCB)[1] == bar(MCB)[2]})
  }
  # B2a
  if(s_==9){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(0,0,0,1/4,0)
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    
    expression_dsc = expression({{bar(DSC)[1] == 0} ~~ plain(and) ~~bar(DSC)[2](italic(k)) == 8 * italic(k)^2 /15})
    expression_mcb <- expression(0 < {bar(MCB)[1] < bar(MCB)[2]})
  }
  # B2b
  if(s_==10){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(0,1/4,0,1/4,0)
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    expression_dsc = expression({{bar(DSC)[1] == 0} ~~ plain(and) ~~bar(DSC)[2](italic(k)) == 8 * italic(k)^2 /15})
    expression_mcb <- expression(0 < {bar(MCB)[1] == bar(MCB)[2]})
  }
  # B3a
  if(s_==11){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(0,0,0,0,k_+1/4)
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    expression_dsc <- expression({bar(DSC)[1](italic(k)) == 2*bar(DSC)[2](italic(k))} == 16 * italic(k)^2 /15)
    expression_mcb <- expression(0 < {bar(MCB)[1] < bar(MCB)[2]})
  }
  # B3b
  if(s_==12){
    beta_=c(rep(0.25,4))
    gamma_=c(0,0,0,k_)
    delta_=c(1/sqrt(15),0,0,0,k_+1/4)
    xi_=c(0,0,k_/2+1/4,0,k_/2+1/4)
    expression_dsc <- expression({bar(DSC)[1](italic(k)) == 2*bar(DSC)[2](italic(k))} == 16 * italic(k)^2 /15)
    expression_mcb <- expression(0 < {bar(MCB)[1] == bar(MCB)[2]})
  }
  
  return(
    list(
      # realization
      beta=beta_, 
      gamma=gamma_,
      # forecasts
      delta=delta_,
      xi=xi_, 
      # meta
      s=s_, k=k_,
      expression_mcb = expression_mcb, 
      expression_dsc = expression_dsc
    )
  )
}

# function to calculate proposition for mean simulation
se_decomp_prop_sim <- function(
    gamma = c(1, 1, 1, 1),
    beta  = rep(0.25, 4),
    parameterforecaster = c(0, 0, 1, 0, 0) # (delta0, delta1, ..., delta4)
){
  
  
  gamma <- gamma %>% t() %>% t()
  varsigma <- 1 / (1 - beta[1]^2)
  
  delta0 <- parameterforecaster[1]
  delta  <- parameterforecaster[2:5]%>% t() %>% t()
  
  dd <- t(delta)%*%(delta)
  dg <- t(delta)%*%(gamma)
  gg <- (t(gamma) %*% gamma)
  
  
  mcb <- delta0^2 + varsigma * ((dd - dg)^2 / dd)
  dsc <- varsigma * (dg^2 / dd)
  unc <- varsigma * gg + 1
  
  list(mcb = as.numeric(mcb), dsc = as.numeric(dsc), unc = as.numeric(unc), varsigma = as.numeric(varsigma))
}

# sets parameters for quantile simulation
parameters_simset_q <- function(s_,k_){
  if(s_==1){
    # para 0<dsc1=dsc2 & mcb1=mcb2 = j
    gamma_ = c(1/4, 1/4, 1/4, 0)
    beta_ = c(rep(.25, 4))
    delta_ = c(0, k_/2+1/4, 0, k_/2+1/4, 0)
    xi_ = c(0, 0, k_/2+1/4, k_/2+1/4, 0)
    
    expression_mcb <- expression({bar(MCB)[1](italic(k)) == bar(MCB)[2](italic(k))})
    expression_dsc = expression({bar(DSC)[1] == bar(DSC)[2]} > 0)
    
    
    
  }
  if(s_==2){
    # para 0<dsc1=dsc2 & 0=mcb1 and mcb2 = j
    gamma_ = c(1/4, 1/4, 1/4, 0)
    beta_ = c(rep(.25, 4))
    delta_ = c(0, 1/4, 0, 1/4, 0)
    xi_ = c(0, 0, k_/2+1/4, k_/2+1/4, 0)
    
    expression_mcb = expression({{bar(MCB)[1](italic(k)) == 0} ~~ plain(and) ~~bar(MCB)[2](italic(k))})     
    expression_dsc = expression({bar(DSC)[1] == bar(DSC)[2]} > 0)
    
  }
  if(s_==3){
    # para dsc1=dsc2 = j & 0 < mcb1 = mcb2 
    gamma_ = c(0, 0, 0, k_)
    beta_ = c(rep(0.25, 4))
    delta_ = c(0, k_/2+1/4, 0, 0, k_/2+1/4) 
    xi_ = c(0, 0, k_/2 + 1/4, 0, k_/2 + 1/4)
    
    expression_mcb = expression({bar(MCB)[1] == bar(MCB)[2]} > 0)
    expression_dsc = expression({bar(DSC)[1](italic(k)) == bar(DSC)[2](italic(k))})
    
    
  }
  if(s_==4){
    # para dsc1=0 and dsc2=j & 0< mcb1=mcb2
    gamma_ = c(0, 0, 0, k_)
    beta_=c(rep(0.25,4))
    delta_=c(1/2, 1/4, 0, 1/4, 0)
    xi_ = c(0, 0, k_/2 + 1/4, 0, k_/2 + 1/4)
    
    # here mcb1 is only approx equal with mcb2 or???
    
    expression_mcb = expression({bar(MCB)[1]  %~~% bar(MCB)[2] > 0})
    expression_dsc = expression({{bar(DSC)[1] == 0} ~~ plain(and) ~~bar(DSC)[2](italic(k))})     
    
  }
  
  return(
    list(
      # realization
      beta=beta_, 
      gamma=gamma_,
      # forecasts
      delta=delta_,
      xi=xi_, 
      # meta
      s=s_, k=k_,
      expression_mcb = expression_mcb, 
      expression_dsc = expression_dsc
    )
  )
}




# function to calculate proposition for quantile simulations

mcb_dsc_checkloss_closed <- function(alpha, beta, gamma, delta, delta0) {
  if (any(alpha <= 0 | alpha >= 1)) stop("alpha must be in (0,1).")
  if (!is.numeric(beta) || length(beta) != 1) stop("beta must be a scalar.")
  if (abs(beta) >= 1) stop("Need |beta| < 1 so varsigma = 1/(1-beta^2) is finite.")
  if (!is.numeric(gamma) || !is.numeric(delta)) stop("gamma and delta must be numeric vectors.")
  if (length(gamma) != length(delta)) stop("gamma and delta must have the same length.")
  if (!is.numeric(delta0) || length(delta0) != 1) stop("delta0 must be a scalar.")
  
  z_alpha <- stats::qnorm(alpha)
  phi     <- stats::dnorm
  Phi     <- stats::pnorm
  
  varsigma <- 1 / (1 - beta^2)
  
  gamma2  <- sum(gamma * gamma)
  sigmaY2 <- varsigma * gamma2 + 1
  sigmaY  <- sqrt(sigmaY2)
  
  delta2 <- sum(delta * delta)
  if (delta2 <= 0) stop("delta^T delta must be > 0 (delta cannot be the zero vector).")
  
  del_gam  <- sum(delta * gamma)
  sigma1X2 <- sigmaY2 - varsigma * (del_gam^2) / delta2
  sigma1X2 <- max(sigma1X2, 0)
  sigma1X  <- sqrt(sigma1X2)
  
  m1 <- -(delta0 + z_alpha)
  
  s12 <- 1 + varsigma * (gamma2 + delta2 - 2 * del_gam)
  if (s12 <= 0) stop("Computed s1^2 <= 0; check your parameters.")
  s1 <- sqrt(s12)
  
  kappa1 <- m1 / s1
  
  MCB1 <- m1 * (Phi(kappa1) + alpha - 1) + s1 * phi(kappa1) - sigma1X * phi(z_alpha)
  DSC1 <- phi(z_alpha) * (sigmaY - sigma1X)
  UNC  <- sigmaY * phi(z_alpha)
  
  list(mcb = MCB1, dsc = DSC1, unc = UNC)
}


# The functions next are needed to find the xi0 for the the 4th q-sim DGP
# Helper: closed-form MCB as a scalar function of delta0
mcb_closed <- function(alpha, param, delta_vec, delta0) {
  mcb_dsc_checkloss_closed(
    alpha  = alpha,
    beta   = unique(param$beta),
    gamma  = param$gamma,
    delta  = delta_vec,
    delta0 = delta0
  )$mcb
}

## Safe solver: find xi0 so MCB1 - MCB2 = 0 
solve_xi0_mcb <- function(alpha, param, delta0_fixed = 0.5,
                          init_width = 1, max_iter = 60) {
  
  # Model 1: delta slopes; Model 2: xi slopes (intercept solved)
  delta_slopes <- param$delta[-1]
  xi_slopes    <- param$xi[-1]
  
  target_mcb1 <- mcb_closed(alpha, param, delta_slopes, delta0_fixed)
  
  g <- function(xi0) {
    target_mcb1 - mcb_closed(alpha, param, xi_slopes, xi0)
  }
  
  # bracket expansion around delta0_fixed
  center <- delta0_fixed
  width  <- init_width
  
  for (i in seq_len(max_iter)) {
    #lo <- center - width
    #hi <- center + width
    
    lo <- max(0, center - width)
    hi <- center + width
    
    glo <- g(lo)
    ghi <- g(hi)
    
    if (is.finite(glo) && is.finite(ghi) && glo * ghi <= 0) {
      return(uniroot(g, interval = c(lo, hi))$root)
    }
    width <- width * 2
  }
  
  # fallback: minimize squared gap if no sign change found
  lo <- center - width
  hi <- center + width
  optimize(function(x) g(x)^2, interval = c(lo, hi))$minimum
}


getxi0 <- function(K, ALPHA, init_width = 1, max_iter = 60) {
  stopifnot(is.numeric(K), length(K) == 1L,
            is.numeric(ALPHA), length(ALPHA) == 1L)
  
  param <- parameters_simset_q(s_ = 4, k_ = K)
  
  solve_xi0_mcb(
    alpha = ALPHA,
    param = param,
    delta0_fixed = param$delta[1],
    init_width = init_width,
    max_iter = max_iter
  )
}
