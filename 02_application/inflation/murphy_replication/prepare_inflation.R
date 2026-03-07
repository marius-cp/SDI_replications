rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("helpers.R")
library(murphydiagram)
data(inflation_mean)  # to check if we get same data as in Ehm et al (2016)


######################################
# Load and prepare Michigan data
######################################

# data: https://data.sca.isr.umich.edu/tables.php  go to nr32, 
# Expected Change in Prices During the Next Year, xls file in last column
# top rows of excel were manually removed 
dat <- readxl::read_xls("redbk32.xls")[-1,]
michigan_dat <- 
  bind_rows(
  dat %>% filter(`Date  of  Survey`=="February") %>% mutate(Q = 1,M=2),
  dat %>% filter(`Date  of  Survey`=="May") %>% mutate(Q = 2,M=5),
  dat %>% filter(`Date  of  Survey`=="August") %>% mutate(Q = 3,M=8),
  dat %>% filter(`Date  of  Survey`=="November") %>% mutate(Q = 4,M=11)
) %>% 
  mutate(
    date0=paste(...2,"Q",Q, sep = ""),
    dt=paste((...2)+1,"Q",Q, sep = ""), 
    stemp_rlz= paste((...2)+1,M,"01",sep = "-") %>% as.Date()
    ) %>% 
  mutate(
    year=...2, 
    quarter = Q
  ) %>% 
  dplyr::select(stemp_rlz,year,quarter,dt,Median)  %>% 
  arrange(stemp_rlz) %>% 
  filter(Median != "NA") %>% 
  mutate(
    value = round(as.numeric(Median),3), 
    id ="michigan") %>% 
  select(stemp_rlz,year,quarter,dt,value, id) %>% 
  filter(stemp_rlz>="1982-08-01")
  
michigan_dat

# check if same as Ehm et al 
all.equal(
  inflation_mean$michigan, 
  michigan_dat %>% filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2014-08-01") %>% pull(value)
)


######################################
# Load and prepare SPF data
######################################

# SPF data: http://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/CPI/
#dat <- readLines("Individual_CPI.csv")

dat <- readxl::read_xlsx("Individual_CPI.xlsx")

spf_dat <- 
dat %>% 
  select(YEAR,QUARTER,CPI3:CPI6) %>% 
  mutate_all(~ ifelse(. == "#N/A", NA, .)) %>% 
  mutate(
    date0=paste(YEAR,"Q",QUARTER, sep = ""),
    dt=paste((YEAR)+1,"Q",QUARTER, sep = ""), 
    M = case_when(
      QUARTER == 1 ~ 2,
      QUARTER == 2 ~ 5,
      QUARTER == 3 ~ 8,
      QUARTER == 4 ~ 11
    ),
    stemp_rlz= paste((YEAR)+1,M,"01",sep = "-")
    ) %>% 
  mutate(
    # Compute mean forecast over next quarters (each forecaster)
    across(c(CPI3, CPI4, CPI5, CPI6), as.numeric),
    mean = c(CPI3+CPI4+CPI5+CPI6) / 4 
    ) %>% 
  filter(mean!="NA") %>% 
  group_by(YEAR,QUARTER) %>% 
  summarise(
    # as done in Ehm et al (2016) and Patton (2020)
    value = median(mean, na.rm = T),
    #  SPF Q4 ahead, we use it to check the robustness 
    # value =  median(CPI6, na.rm = T) # results will remain robust 
    ) %>% 
  mutate(
    dt=paste((YEAR)+1,"Q",QUARTER, sep = ""), 
    M = case_when(
      QUARTER == 1 ~ 2,
      QUARTER == 2 ~ 5,
      QUARTER == 3 ~ 8,
      QUARTER == 4 ~ 11
    ),
    stemp_rlz= paste((YEAR)+1,M,"01",sep = "-") %>% as.Date(), 
    year=YEAR,quarter=QUARTER, 
    id ="spf"
  ) %>% ungroup() %>% 
  select(stemp_rlz,year,quarter,dt, value, id)

spf_dat

# allmost eual !!!
bind_cols(
  spf_dat   %>% 
   filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2014-08-01"),# %>%
    #pull(spf),# %>% length()
 inflation_mean #%>% pull(spf)# %>% length()
) %>% mutate(dif = value - spf) %>% arrange(dif) %>% 
  select(year, quarter,dif)


all.equal(
  spf_dat   %>% 
    filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2014-08-01") %>%
    pull(value),# %>% length()
  inflation_mean %>% pull(spf)# %>% length()
)


######################################
# Load and prepare CPI data
######################################
# data: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/cpi
dat <- readxl::read_xlsx("cpiQvMd.xlsx")
dat <- dat %>% select(DATE,CPI23Q3) %>% data.frame() 
dat


# some cleaning
dat[, 2] <- gsub(",", ".", dat[,2])
dat[, 2][dat[, 2] == ""] <- NA
dat[,2] <- as.numeric(dat[,2])
dat <- dat[!is.na(dat[,2]), ]
names(dat)[2] <- "val"
   
dat$y <- substr(dat[,1], 1, 4)
dat$m <- as.numeric(substr(dat[,1], 6, 7))
dat$q <- (dat$m %in% 1:3) + 2*(dat$m %in% 4:6) + 3*(dat$m %in% 7:9) + 4*(dat$m %in% 10:12)
dat$dt <- paste0(dat$y, "Q", dat$q)
dat$stemp_rlz <-  as.Date(paste(dat$y,dat$m,01,sep="-"))
dat <- subset(dat, select = c("dt", "val","y","m","q"))
dat %>% group_by(y,q) %>% summarise(cpi=mean(val)) %>% pull(cpi)

all.equal(
aggregate(dat$val, by = list(dat$dt), mean) %>% pull(x), 
dat %>% group_by(y,q) %>% summarise(cpi=mean(val)) %>% pull(cpi)
)

cpi_dat <- 
  dat %>% group_by(y,q) %>% 
  summarise(cpi=mean(val), m = median(m)) %>% 
  ungroup() %>% # important !!!!
  mutate(
    log_cpi = log(cpi),
    value = c(rep(NA,4),100*diff(log_cpi, lag=4)),
    dt=paste((y),"Q",q, sep = ""), 
    stemp_rlz= paste((y),m,"01",sep = "-") %>% as.Date(),
    year = as.double(y)-1, quarter=q,id="rlz"
  ) %>% 
  select(stemp_rlz,year, quarter, dt, value, id)

cpi_dat

all.equal(
  cpi_dat %>% filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2014-08-01") %>% pull(value),
  inflation_mean %>% pull(rlz)
  )

check <-
bind_cols(
  cpi_dat %>% filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2014-08-01") ,
  inflation_mean 
) %>% 
  mutate(dif=value-rlz)

######################################
# Joint data frame
######################################


inflation_mean_new <- 
bind_rows(
  cpi_dat %>% filter(stemp_rlz >= "1982-08-01"& stemp_rlz <= "2023-05-01"),
  spf_dat %>% filter(stemp_rlz >= "1982-08-01"& stemp_rlz <= "2023-05-01"),
  michigan_dat %>% filter(stemp_rlz >= "1982-08-01" & stemp_rlz <= "2023-05-01")
) %>% ungroup() %>% arrange(stemp_rlz,id) %>% 
  pivot_wider(
  names_from = id, 
  values_from = c(value)
) %>% 
  select(!year) %>% 
  select(stemp_rlz,quarter,dt,spf,michigan,rlz)

inflation_mean_new

saveRDS(inflation_mean_new, "inflation_mean.rds")
