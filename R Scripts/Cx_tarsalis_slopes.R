## libraries----------------------------------------------------------------
library(wsyn)

# libraries for data viz
library(ggplot2)
library(patchwork)


## load and prep Cx tarsalis timeseries data--------------------------------

times<-1:45

############# cx tarsalis data ############# 

dat<-as.matrix(read.csv("neon_cx_tarsalis.csv", header = T))
dat2<-dat[,-1] ## removes tBin
dat2<-dat2[,-3] ## removes JORN
 
dat2<-t(dat2) # transpose data

# clean data:  
x<-cleandat(dat2,times,3)$cdat
# 3: time series are (individually) linearly detrended and de-meaned, and variances are standardized to 1


## Cx wmf-------------------------------------------------------------------
moswmf<-wmf(x,times = times)

# get value matrix
gwmf<-abs(moswmf$value)
# gwmf


## Cx wmf plot syn ---------------------------------------------------------
plotmag(moswmf)




## Cx wmf extract timescales------------------------------------------------
# 12 month time scale
gwmf12<-gwmf[12,]
gwmf12<-gwmf12[-46]

gwmf12

# 8-14 month timescales average
gwmf8_14<-colMeans(gwmf[8:14,], na.rm = T)
gwmf8_14<-gwmf8_14[-46]

gwmf8_14


## Cx wmf slopes------------------------------------------------------------
# 12 month time scale
out1<-lm(gwmf12 ~ times)
slope_real<- out1$coefficients[2]

# 8-14 month timescales average
out2<-lm(gwmf8_14 ~ times)
slope_real2 <- out2$coefficients[2]

slopes_wmf <- cbind(slope_real, slope_real2)


## Surr data generation-----------------------------------------------------
## Get synchrony preserving aaft surrogate data of your mosquito data
nsurrogs <- 1 # for code check
surrtype <- "aaft"
syncpres <- TRUE
sur <- surrog(x,nsurrogs,surrtype,syncpres)
## this is a demo of the used code, we generated 10,000 surr data sets using a looped setup

# with x being the prep. mosquito data:

    # dat<-as.matrix(read.csv("neon_cx_tarsalis.csv", header = T))
    # dat2<-dat[,-1]
    # dat2<-dat2[,-3]
    # dat2<-t(dat2)
    # x<-cleandat(dat2,times,3)$cdat




## Get sur wmf slopes-------------------------------------------------------

df_slope_wmf <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(df_slope_wmf) <- c("slope_surrog","slope_surrog2")

for (i in sur) {
  
  # surrogate data
  x <- i 
  
  # times<-1:45

  moswmf_s<-wmf(x,times = times)
  
  gwmfs<-abs(moswmf_s$value)
  
  # get timeseries values:
  
  gwmfs12<-gwmfs[12,] 
  gwmfs12<-gwmfs12[-46] 
  
  gwmfs8_14<-colMeans(gwmfs[8:14,], na.rm = T)
  gwmfs8_14<-gwmfs8_14[-46]
  
  # get slopes:
  
  out1s<-lm(gwmfs12 ~ times)
  slope_surrog <- out1s$coefficients[2]
  
  out2s<-lm(gwmfs8_14 ~ times)
  slope_surrog2 <- out2s$coefficients[2]
  
  slopes2 <- cbind(slope_surrog, slope_surrog2)
  df_slope_wmf <- rbind(df_slope_wmf, slopes2)
  
}



## -------------------------------------------------------------------------
## read slopes for all 10000 surrogates
slp_Swmf <- read.csv("output/surr_slopes/wmf_Surrogate_slopes_aaft_10000_140524.csv", stringsAsFactors = F)
rslp_12_wmf <- slopes_wmf[1,2]

h1a <- ggplot(slp_Swmf, aes(slope_surrog))+
  geom_histogram(fill="lavenderblush4", alpha=.8) +
  geom_vline(xintercept=0, linewidth = 0.5, linetype = 2, color = "black") + 
  geom_vline(xintercept=slopes_wmf[1,1], linewidth = 1, color = "orangered3") +
  labs(title = 'wmf()', 
       subtitle="12 month timescale") +
  theme_bw() 

h1b <- ggplot(slp_Swmf, aes(slope_surrog2))+
  geom_histogram(fill="lavenderblush4", alpha=.8) +
  geom_vline(xintercept=0, linewidth = 0.5, linetype = 2, color = "black") + 
  geom_vline(xintercept=slopes_wmf[1,2], linewidth = 1, color = "orangered3") +
  labs(subtitle="8-14 month avg. timescale") +
  theme_bw() 

h1a + h1b


## plot side by side for comparison --------------------------------------------
par(mar = c(4, 4, .1, .1))
plotmag(moswmf_s)
plotmag(moswmf)




## Cx wpmf------------------------------------------------------------------

moswpmf<-wpmf(x,times = times,sigmethod = "aaft")

# get value matrix
gwpmf<-abs(moswpmf$value)
gwpmf

plotmag(moswpmf)


## Cx wpmf extract timescales-----------------------------------------------
# 12 month time scale
gwpmf12<-gwpmf[12,]
gwpmf12<-gwpmf12[-46]

gwpmf12

# 8-14 month timescales average
gwpmf8_14<-colMeans(gwpmf[8:14,], na.rm = T)
gwpmf8_14<-gwpmf8_14[-46]

gwpmf8_14


## Cx wpmf slopes-----------------------------------------------------------
# 12 month time scale
out1<-lm(gwpmf12 ~ times)
slope_real<- out1$coefficients[2]

# 8-14 month timescales average
out2<-lm(gwpmf8_14 ~ times)
slope_real2 <- out2$coefficients[2]

slopes_wpmf <- cbind(slope_real, slope_real2)


## Get sur wpmf slopes------------------------------------------------------

df_slope_wpmf <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(df_slope_wpmf) <- c("slope_surrog","slope_surrog2")

for (i in sur) {
  
  # surrogate data
  x <- i 
  
  # times<-1:45

  moswpmf_s<-wpmf(x,times = times,sigmethod = "aaft")
  
  gwpmfs<-abs(moswpmf_s$value)
  
  # get timeseries values:
  
  gwpmfs12<-gwpmfs[12,] 
  gwpmfs12<-gwpmfs12[-46] 
  
  gwpmfs8_14<-colMeans(gwpmfs[8:14,], na.rm = T)
  gwpmfs8_14<-gwpmfs8_14[-46]
  
  # get slopes:
  
  out1s<-lm(gwpmfs12 ~ times)
  slope_surrog <- out1s$coefficients[2]
  
  out2s<-lm(gwpmfs8_14 ~ times)
  slope_surrog2 <- out2s$coefficients[2]
  
  slopes2 <- cbind(slope_surrog, slope_surrog2)
  df_slope_wpmf <- rbind(df_slope_wpmf, slopes2)
  
}



## -------------------------------------------------------------------------
## read in slopes for 10000 surrogate data sets
slp_Swpmf <- read.csv("output/surr_slopes/Surr_slopes_aaft_10k.csv", stringsAsFactors = F)

h2a <- ggplot(slp_Swpmf, aes(surSlp12))+
  geom_histogram(fill="lightsteelblue4", alpha=.8) +
  geom_vline(xintercept=0, linewidth = 0.5, linetype = 2, color = "black") + 
  geom_vline(xintercept=slopes_wpmf[1,1], linewidth = 1, color = "blue4") +
  labs(title = 'wpmf()', 
       subtitle="12 month timescale") +
  theme_bw() 

h2b <- ggplot(slp_Swpmf, aes(surSlp8_14))+
  geom_histogram(fill="lightsteelblue4", alpha=.8) +
  geom_vline(xintercept=0, linewidth = 0.5, linetype = 2, color = "black") + 
  geom_vline(xintercept=slopes_wpmf[1,2], linewidth = 1, color = "blue4") +
  labs(subtitle="8-14 month avg. timescale") +
  theme_bw() 

h2a + h2b


##  plot side by side for for comparison-------------------------------------------
par(mar = c(4, 4, .1, .1))
plotmag(moswpmf_s)
plotmag(moswpmf)


