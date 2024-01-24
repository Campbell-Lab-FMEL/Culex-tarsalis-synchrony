### wavelet linear model 

# load packages
library(wsyn)
library(stringr)
library(tidyverse)

## libraries for data viz
# library(ggplot2)
# library(patchwork)


# define workspace: 
setwd("C:/Users/amely/OneDrive - University of Florida/B_Campbell_Lab/2023_Synchrony_CxTarsalis/2024_Cx_syn")

workDir <- file.path(".")
inDir <- file.path(workDir, "input")
outDir <- file.path(workDir, "output")
if(!dir.exists(outDir)) dir.create(outDir)

#########################################################

# check content of data folder
list.files(inDir)

# define environmental variables included in analyses
env <- c("prec", "tmin", "tmax", "rmin", "rmax") 


#########################################################
### Prepare data for analyses

dat.list <- list() # store all data in this list

# pattern in env. variable data file names
p <- "_for_wsyn_times_45.csv"

## list all data files in data folder:
input.files <- list.files(inDir)

## loop through all files
# remove unnecessary columns, transpose & save as matrix objects in the dat.list list

for(f in input.files){
  
  # load data
  d <- read.csv(file.path(inDir, f))
  
  # if file = species data
  if(str_detect(f, "cx_tar")){
    
    v.name <- "cx_tar"

    d1 <- d %>% 
      arrange(tBin) %>% 
      select(-c(tBin)) %>% 
      t() %>% 
      as.matrix()
    
    dat.list[[v.name]] <- d1
    
  }
  
  # if file = env data
  if(str_detect(f, p)){
    
    v.name <- f %>% str_remove(p) %>% word(-1, sep = "_")

    d1 <- d %>% 
      arrange(tBin) %>% 
      select(-c(year, month, tBin)) %>% 
      t() %>% 
      as.matrix()
    
    dat.list[[v.name]] <- d1
    
  }
  
  # save data as .RDS object
  if(f == last(input.files)){
    
    # rearrange list objects:
    dat <- dat.list[c("cx_tar", env)] 
    saveRDS(dat, file.path(outDir, "cxtar_wsyn45_data.RDS"))
    
    }
  
}

#########################################################
### prep for wavelet linear models (wlm)

## wlm()
# define timesteps
times <- 1:45
# clean data
dat <- lapply(FUN=function(x){cleandat(x,times,3)$cdat},X=dat)

# define wlm structure
resp <- 1
pred <- 2:6 # will loop through these 
norm <- "powall"

## wlmtest()
# define nrand
n <- 10000
sigmethod <- "aaft"

## bandtest() 
# timescale of interest
blong<-c(11,13) # ~annual timescale, 11-13 months

res.list <- list()

### ### ### ### ### ### ### ### ###
### define wlm & test for significant relationship between response and predictor

# fit and test wlm for each environmental variable:
# res represents the significance of the variable dropped (!!!)

for(p in pred){
  # control:
  cat("env. predictor:", p, "\n", sep = " ")
  
  # define wlmobj
  dat.wlm <- wlm(dat, times, resp, p, norm)
  (print)
  
  # wlmtest
  res <- wlmtest(dat.wlm, drop = 2, sigmethod, nrand= n) # drop dat.wlm$dat[2], i.e. the env covariate
  res <- bandtest(res, blong)

  k <- env[p-1]

  # saves plotrank plots in the out folder
  plotrank(res)
  plotrank(res,
           filename = file.path(outDir, str_c("cxtar_wlmtest_", k, "_", n, sep = "")))

  res.list[[k]] <- res
    if(p == last(pred)){saveRDS(res.list, file.path(outDir, "cxtar_wsyn45_wlmtest.RDS"))}

}

