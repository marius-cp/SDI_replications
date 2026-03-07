rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(lubridate)
library(tidyverse)
library(midasr)
library(rugarch)
library(doParallel)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# load cleaned emini data
dat <- readRDS("data/emini_clean.rds") %>% mutate(Loc_date=stday)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# construct HAR predictors
df <- dat %>%
  mutate(
    y = sqrt(RV_5min),
    Lag1 = lag(y, 1),
    Lag5 = lag(zoo::rollapply(y, 5, mean, align = "right", fill = NA), 1),
    Lag20 = lag(zoo::rollapply(y, 22, mean, align = "right", fill = NA), 1)
  ) %>%
  ungroup()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Split into in-sample and out-of-sample
is <- df %>% filter(Loc_date > "2000-01-01", Loc_date < "2008-01-01")
os <- df %>% filter(Loc_date > "2008-01-01", Loc_date < "2022-01-01")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# GARCH forecasts
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "std"
  )
fit_garch <- ugarchfit(spec = spec_garch, data = is$ret_oc)
spec_garch_fixed <- rugarch:::.setfixed(spec_garch, as.list(coef(fit_garch)))
filtered_garch_os <- ugarchfilter(spec = spec_garch_fixed, data = os$ret_oc)
f <- function(x, fit) fGarch::qstd(p=x, mean=0, sd=1, nu=fit@fit$coef["shape"]) 
mu_t    <- fitted(filtered_garch_os)                
sigma_t <- matrix(sigma(filtered_garch_os), ncol=1) 
VaR_sgarch_t <- mu_t + sigma_t %*% sapply(0.01, f, fit_garch)
VaR_sgarch_t_5 <- mu_t + sigma_t %*% sapply(0.05, f, fit_garch)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# GJR forecasts
spec_gjr <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "std"
)
fit_gjr <- ugarchfit(spec = spec_gjr, data = is$ret_oc)
spec_gjr_fixed <- rugarch:::.setfixed(spec_gjr, as.list(coef(fit_gjr)))
filtered_gjr_os <- ugarchfilter(spec = spec_gjr_fixed, data = os$ret_oc)
gjr_mu_t    <- fitted(filtered_gjr_os)                
gjr_sigma_t <- matrix(sigma(filtered_gjr_os), ncol=1) 
VaR_gjrgarch_t <- gjr_mu_t + gjr_sigma_t %*% sapply(0.01, f, fit_gjr)
VaR_gjrgarch_t_5 <- gjr_mu_t + gjr_sigma_t %*% sapply(0.05, f, fit_gjr)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# HAR forecasts
fit_har <- lm(y ~ Lag1 + Lag5 + Lag20, data = is)
pred_har <- predict(fit_har, newdata = os)
VaR_HAR <- as.numeric(pred_har) * qnorm(0.01)
VaR_HAR_5 <- as.numeric(pred_har) * qnorm(0.05)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# MIDAS forecasts
midasdat <- bind_rows(is, os)
midas_model <- midas_r(
  rv ~ mls(rv, 1:20, 1, nealmon),
  data = list(rv = sqrt(is$RV_5min)),
  start = list(rv = c(0, 0, 0)),
  weight_gradients = list()
)
midasfcast <- average_forecast(
  list(midas_model),
  data = list(rv = sqrt(midasdat$RV_5min)),
  insample = 1:nrow(is),
  outsample = (nrow(is) + 1):nrow(midasdat),
  type = "fixed",
  show_progress = FALSE
)
VaR_MIDAS <-  midasfcast$forecast * qnorm(0.01)
VaR_MIDAS_5 <- midasfcast$forecast * qnorm(0.05)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# HS  forecasts (only for VaR)
window<-250
HSq <- data.frame(date=df$Loc_date, HSfcast=rep(NA_real_,nrow(df)))   
for (t in (window+1):nrow(df)) {
  # use the previous 250 daily returns: (t-250) ... (t-1)
  window_ret <- df$ret_oc[(t - window):(t - 1)]
  # alpha-quantile of past returns = 1-step-ahead HS VaR for day t
  HSq$HSfcast[t] <- stats::quantile(window_ret, probs = 0.01, na.rm = TRUE, type=1)
}
HSq <- HSq %>% filter(date > "2008-01-01", date < "2022-01-01")

HSq_5 <- data.frame(date = df$Loc_date, HSfcast = rep(NA_real_, nrow(df)))   
for (t in (window + 1):nrow(df)) {
  window_ret <- df$ret_oc[(t - window):(t - 1)]
  HSq_5$HSfcast[t] <- stats::quantile(
    window_ret, 
    probs = 0.05, 
    na.rm = TRUE, 
    type = 1
  )
}
HSq_5 <- HSq_5 %>% filter(date > "2008-01-01", date < "2022-01-01")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Collect vola forecasts in one tibble
fcasts_vola <- tibble(
  Symbol = "Emini",
  date = os$Loc_date,
  # realized variance:
  rv5 = os$RV_5min,
  # open-close returns:
  ret = os$ret_oc,
  # forecast of realized variance
  HAR = as.numeric(pred_har)^2,
  MIDAS = as.numeric(midasfcast$forecast)^2,
  GJR = as.numeric(sigma(filtered_gjr_os))^2,
  GARCH = as.numeric(sigma(filtered_garch_os))^2
  ) 
saveRDS(fcasts_vola, "data/emini_fcasts_vola.rds")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Collect 1% VaR forecasts in one tibble
fcasts_var_1 <- tibble(
  Symbol = "Emini",
  date = os$Loc_date,
  # realized variance:
  rv5 = os$RV_5min,
  # open-close returns:
  ret = os$ret_oc,
  # forecast of realized variance
  HAR = as.numeric(VaR_HAR),
  MIDAS = as.numeric(VaR_MIDAS),
  GJR = as.numeric(VaR_gjrgarch_t),
  GARCH = as.numeric(VaR_sgarch_t), 
  HS = as.numeric(HSq$HSfcast)
) 
saveRDS(fcasts_var_1, "data/emini_fcasts_VaR_1.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Collect 5% VaR forecasts in one tibble
fcasts_var_5 <- tibble(
  Symbol = "Emini",
  date   = os$Loc_date,
  rv5    = os$RV_5min,
  ret    = os$ret_oc,
  HAR    = as.numeric(VaR_HAR_5),
  MIDAS  = as.numeric(VaR_MIDAS_5),
  GJR    = as.numeric(VaR_gjrgarch_t_5),
  GARCH = as.numeric(VaR_sgarch_t_5),
  HS     = as.numeric(HSq_5$HSfcast)
)
saveRDS(fcasts_var_5, "data/emini_fcasts_VaR_5.rds")
