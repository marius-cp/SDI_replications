rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(lubridate)
library(tidyverse)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# emini raw ----
emini_raw <- readRDS("./data/index_futures_ESc1_1min_2000_2022.rds")
min5 <- seq(0,60,5)
min10 <- seq(0,60,10)
min30 <- seq(0,60,30)
sc <- 100

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cleaning step 1 ----
# extract infos 
emini_clean1 <- 
  emini_raw %>%   
  mutate(
    GMT_dttm = with_tz(`Date-Time`, tzone ="GMT"),
    Loc_dttm = GMT_dttm+hours(`GMT Offset`),
    Loc_date = date(Loc_dttm),
    Loc_time = format(Loc_dttm, format = "%H:%M:%S"), 
    tradingday = Loc_dttm-hours(17), 
    tday_count = with(rle(as.numeric(date(tradingday))), rep(seq_along(lengths), times = lengths)), # counts trading days 
    tday_time = format(tradingday, format = "%H:%M:%S"), # hms of the trading day 
    weekday  =  weekdays(Loc_date),
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cleaning step 2----
# trading hours 
emini_clean2 <- 
  emini_clean1 %>% 
  dplyr::filter(
    weekday != "Saturday" ,
    !(Loc_time >= "15:15:00" & Loc_time < "15:30:00"),  # maintenance
    !(Loc_time >= "16:00:00" & Loc_time < "17:00:00") , # stop trading
    !(weekday == "Friday" & Loc_time > "16:00:00"), 
    !(weekday == "Sunday" & Loc_time < "17:00:00")
  ) %>% 
  mutate(
    t.price = Last, 
    t.price.was.na =ifelse(is.na(t.price),1,0)
  ) %>% 
  dplyr::select(Loc_dttm, Loc_date, tday_count, tday_time, Volume, t.price, t.price.was.na, Last)

nacheck <- 
  emini_clean2 %>% 
  group_by(tday_count) %>% 
  summarise(
    Loc_date=min(Loc_date),
    nas_per_day=sum(is.na(t.price)), 
    priceobs_per_day=sum(!is.na(t.price)),
    total = nas_per_day+priceobs_per_day
  ) %>% 
  mutate(badday = as.numeric(priceobs_per_day < 100)) %>% 
  select(tday_count,badday)

nacheck$badday %>% mean()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cleaning step 3 ----
# filling up NA prices
emini_clean3 <- 
  left_join(
  emini_clean2,nacheck
  ) %>% 
  filter(
    # get rid of bad days 
    badday == 0
  ) %>% 
  fill(t.price, .direction = "down")   # carry last observed value forward
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cleaning step 4 ----
# calculate m-minute wise returns
emini_clean4 <- 
  emini_clean3 %>% 
  group_by(tday_count) %>% 
  mutate(
    ret_1min =  c(NA,sc*diff(log(t.price))), 
    ret_5min = ifelse(lubridate::minute(lubridate::hms(tday_time))%in%min5,sc*log(t.price/dplyr::lag(t.price, 5)),NA ), # when sum is wanted: r_{t,t-5}=r_{t-4}+...+r_t
    ret_10min = ifelse(lubridate::minute(lubridate::hms(tday_time))%in%min10,sc*log(t.price/dplyr::lag(t.price, 10)),NA ),
    ret_30min = ifelse(lubridate::minute(lubridate::hms(tday_time))%in%min30,sc*log(t.price/dplyr::lag(t.price, 30)),NA )
  ) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# cleaning step 5 ----
#construct measures
emini_clean5 <- 
  emini_clean4 %>% 
  group_by(tday_count) %>% 
  summarise(
    stday = min(Loc_dttm),
    tot_t.prices = sum(!is.na((t.price))),
    P.open = t.price[which(tday_time == min(tday_time))], # first price 
    P.close = t.price[which(tday_time == max(tday_time))], # last price 
    RV_1min = (sum(ret_1min^2, na.rm = T)), # rel variance 1min 
    RV_5min = (sum(ret_5min^2, na.rm = T)),  # rel variance 5min
    RV_10min = (sum(ret_10min^2, na.rm = T)),  # rel variance 10min
    RV_30min = (sum(ret_30min^2, na.rm = T)),  # rel variance 10min
    RVplus_5min = sum(ret_5min^2 * as.numeric(ret_5min>0), na.rm = T), 
    RVminus_5min = sum(ret_5min^2 * as.numeric(ret_5min<0), na.rm = T)
  )%>% 
  ungroup() %>% 
  mutate(
    ret_oc =  sc*log(P.close/P.open), # open-close return
    ret_cc = c(NA, sc*diff(log(P.close))) # close-close return 
  ) %>% 
  filter(
    # drop first day manually due to NA at the beginning
    stday>"2000-01-03"
    )


saveRDS(emini_clean5,"data/emini_clean.rds")

