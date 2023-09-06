setwd("/blue/clord/amelybauer/Cx_tarsalis_syn")

library(wsyn)
library(stringr)
library(doParallel)
#library(purrr)

####### Preparations #############################################################

num_cores <- 28  # Adjust as needed - set to 2 when running on laptop!

## get time for output name:
hms <- format(Sys.time(), format="%H%M%S")

## seq length of timeseries:
times<-1:45

out_path <- str_c("data/output/", hms, "/")
if (!dir.exists(out_path)) dir.create(out_path)

####### Mosquito/Surrogate data: ###############

## prep Cx tarsalis data:
  dat<-as.matrix(read.csv("data/neon_cx_tarsalis.csv", header = T))
  dat2<-dat[,-1] ## removes tBin
  dat2<-dat2[,-3] ## removes JORN
  
  dat2<-t(dat2) # transpose data
  
  # clean data:  
  # 3: time series are (individually) linearly detrended and de-meaned, and variances are standardized to 1
x<-cleandat(dat2,times,3)$cdat
  

## Get surrogate data
# synchrony preserving aaft surrogates of mosquito data

   nsurrogs <- 10000
   surrtype <- "aaft"
   syncpres <- TRUE

  sur <- surrog(x,nsurrogs,surrtype,syncpres)

# save generated surrogate data to be able to check them later: 
sur.data.name <- str_c(out_path, "Cxtarsalis_surrs_", nsurrogs, "_",  surrtype, "_syncpres.csv", stringsAsFactors = F)
write.csv(sur, sur.data.name, row.names=F)
  
  
############# Cx tarsalis wmf/wpmf ############################################### 

### get wmf object and extract timescales ######
moswmf<-wmf(x,times = times)

## extract timescales
  # get value matrix
  gwmf<-abs(moswmf$value)
  
  # 12 month time scale
  gwmf12<-gwmf[12,]
  gwmf12<-gwmf12[-46]
  
  # 8-14 month timescales average
  gwmf8_14<-colMeans(gwmf[8:14,], na.rm = T)
  gwmf8_14<-gwmf8_14[-46]

## calculate slopes:
  # 12 month time scale
  out1<-lm(gwmf12 ~ times)
  slope_real<- out1$coefficients[2]
  
  # 8-14 month timescales average
  out2<-lm(gwmf8_14 ~ times)
  slope_real2 <- out2$coefficients[2]
  
real_slopes_wmf <- cbind(slope_real, slope_real2)

file.name.Cx1 <- str_c(out_path, "Slopes_wmf_Cxtarsalis_", hms, ".csv")
write.csv(real_slopes_wmf, file = file.name.Cx1)

### get wpmf object and extract timescales #####
moswpmf<-wpmf(x,times = times,sigmethod = "aaft")

# get value matrix
gwpmf<-abs(moswpmf$value)

## extract timescales
  # 12 month time scale
  gwpmf12<-gwpmf[12,]
  gwpmf12<-gwpmf12[-46]
  
  # 8-14 month timescales average
  gwpmf8_14<-colMeans(gwpmf[8:14,], na.rm = T)
  gwpmf8_14<-gwpmf8_14[-46]

## calculate slopes:
  # 12 month time scale
  out1<-lm(gwpmf12 ~ times)
  slope_real<- out1$coefficients[2]
  
  # 8-14 month timescales average
  out2<-lm(gwpmf8_14 ~ times)
  slope_real2 <- out2$coefficients[2]
  
real_slopes_wpmf <- cbind(slope_real, slope_real2)

file.name.Cx2 <- str_c(out_path, "Slopes_wpmf_Cxtarsalis_", hms, ".csv")
write.csv(real_slopes_wpmf, file = file.name.Cx2)



############# Surrogate wmf/wpmf ################################################# 

####### Surrogate wmfs #############################

#Set the number of cores for parallel processing
#num_cores <- detectCores(logical = FALSE)

# Create a parallel cluster and register it
cl <- makeCluster(num_cores)
#cl <- makeCluster(floor(0.75 * num_cores)
registerDoParallel(cl)


# function for surrogate wmfs
  # create wmf object, get timescales and calculate slopes
  # returns list of slopes and wmf objects 

wmf_iteration <- function(i) {
  
  xs <- i  # Surrogate data
  
  moswmf_s<-wmf(xs,times = times)
  gwmfs<-abs(moswmf_s$value)
  
  # get timeseries values:
  gwmfs12<-gwmfs[12,] 
  gwmfs12<-gwmfs12[-46] 
  
  gwmfs8_14<-colMeans(gwmfs[8:14,], na.rm = T)
  gwmfs8_14<-gwmfs8_14[-46]
  
  # get slopes:
  
  out1s<-lm(gwmfs12 ~ times)
  slope_surrog <- out1s$coefficients[2]
  names(slope_surrog) <- "slope_surrog"
  
  out2s<-lm(gwmfs8_14 ~ times)
  slope_surrog2 <- out2s$coefficients[2]
  names(slope_surrog2) <- "slope_surrog2"
  
  df_slope <- list(slope_surrog = slope_surrog, slope_surrog2 = slope_surrog2, gwmfs_list = gwmfs)
  
  # Return df_slope and gwmfs as separate lists
  return(list(df_slope))
}


# Apply the function in parallel
results <- foreach(j = 1:length(sur), .combine = c, .multicombine = TRUE, .init = list(), .packages = "wsyn") %dopar% {
                     wmf_iteration(sur[[j]])
                   }

# Stop the parallel cluster
stopCluster(cl)


# Extract the slopes dataframe and the list of gwpmfs from the results
# Combine the slope results into a data frame
sur1 <- do.call(rbind, lapply(results, `[[`, 1))
sur2 <- do.call(rbind, lapply(results, `[[`, 2))
sur_slopes_wmf <- data.frame(sur1, sur2)

# write list of all value matrices of surrogate wmfs
sur_wmf_list  <- lapply(results, `[[`, 3)


# Save the resulting sur slope data frame
file.name.wmf <- str_c(out_path, "Slopes_wmf_sur_", surrtype, "_", length(sur), "_", hms, ".csv")
write.csv(sur_slopes_wmf, file = file.name.wmf)


####### Surrogate wpmfs #############################


#Set the number of cores for parallel processing
#num_cores <- detectCores(logical = FALSE)

# Create a parallel cluster and register it
cl <- makeCluster(num_cores)
#cl <- makeCluster(floor(0.75 * num_cores)
registerDoParallel(cl)


# function for surrogate wpmfs
  # create wpmf object, get timescales and calculate slopes
  # returns list of slopes and wpmf objects

wpmf_iteration <- function(i) {

  xs <- i  # Surrogate data

  moswpmf_s<-wpmf(xs,times = times,sigmethod = "aaft")
  gwpmfs<-abs(moswpmf_s$value)

  # get timeseries values:
  gwpmfs12<-gwpmfs[12,]
  gwpmfs12<-gwpmfs12[-46]

  gwpmfs8_14<-colMeans(gwpmfs[8:14,], na.rm = T)
  gwpmfs8_14<-gwpmfs8_14[-46]

  # get slopes:

  out1s<-lm(gwpmfs12 ~ times)
  slope_surrog <- out1s$coefficients[2]
  names(slope_surrog) <- "slope_surrog"

  out2s<-lm(gwpmfs8_14 ~ times)
  slope_surrog2 <- out2s$coefficients[2]
  names(slope_surrog2) <- "slope_surrog2"

df_slope <- list(slope_surrog = slope_surrog, slope_surrog2 = slope_surrog2, gwpmfs_list = gwmfs)

# Return df_slope and gwmfs as separate lists
return(list(df_slope))
}

# Apply the function in parallel
results <- foreach(j = 1:length(sur), .combine = c, .multicombine = TRUE, .init = list(), .packages = "wsyn") %dopar% {
  wmf_iteration(sur[[j]])
}

# Stop the parallel cluster
stopCluster(cl)


# Extract the slopes dataframe and the list of gwpmfs from the results
# Combine the slope results into a data frame
sur1 <- do.call(rbind, lapply(results, `[[`, 1))
sur2 <- do.call(rbind, lapply(results, `[[`, 2))
sur_slopes_wpmf <- data.frame(sur1, sur2)

# write list of all value matrices of surrogate wmfs
sur_wpmf_list  <- lapply(results, `[[`, 3)


# Save the resulting sur slope data frame
file.name.wpmf <- str_c(out_path, "Slopes_wpmf_sur_", surrtype, "_", length(sur), "_", hms, ".csv")
write.csv(sur_slopes_wpmf, file = file.name.wpmf)



############# save objects ############################################### 

save_for_inspection <- c("gwmf", "gwpmf", "real_slopes_wmf", "real_slopes_wpmf", "sur", "sur_slopes_wmf", "sur_wmf_list", "sur_slopes_wpmf", "sur_wpmf_list")

#save_for_inspection <- c("gwmf", "gwpmf", "real_slopes_wmf", "real_slopes_wpmf", "sur_slopes_wmf", "sur_wmf_list")

r_name <- str_c(out_path, "Cx_sur_", hms, ".RData")
save(list = save_for_inspection, 
     file = print(r_name))

