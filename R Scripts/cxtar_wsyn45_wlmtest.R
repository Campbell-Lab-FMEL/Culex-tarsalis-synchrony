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
### Set up wavelet linear models (wlm)

# define timesteps
times <- 1:45
# clean data
dat <- lapply(FUN=function(x){cleandat(x,times,3)$cdat},X=dat)

# generate wlm object:
resp <- 1
pred <- 2:6  
norm <- "powall"

dat.wlm <- wlm(dat, times, resp, pred, norm)


##################################################################################
### compare wlm - test for significant relationship between response and predictor

# timescale of interest
blong<-c(11,13) # ~annual timescale, 11-13 months

res.list <- list()

# wlmtest for each environmental variable: 
for(i in pred){
  
  drop <- pred[!pred %in% i] 
  
  # control:
  cat("env. predictor:", i, "\n", sep = " ")
  cat("drop:", drop, "\n", sep = " ")
  
  sigmethod<-"aaft"
  n = 100 # should at least be 1000, better 10000
  
  res <- wlmtest(dat.wlm, drop, sigmethod, nrand= n)
  res <- bandtest(res, blong)
  
  k <- env[i-1]
  
  # saves plotrank plots in the out folder
  plotrank(res, 
           filename = file.path(outDir, str_c("cxtar_wlmtest_", k, "_", n, sep = "")))
  
  
  res.list[[k]] <- res
    if(i == last(pred)){saveRDS(res.list, file.path(outDir, "cxtar_wsyn45_wlmtest.RDS"))}
  
}

