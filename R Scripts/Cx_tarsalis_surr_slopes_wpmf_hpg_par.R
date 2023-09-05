setwd("/blue/clord/amelybauer/Cx_tarsalis_syn")

library(wsyn)
library(stringr)
library(doParallel)


## get time for output end:
hms <- format(Sys.time(), format="%H%M%S")

############# Prep Cx tarsalis data ############# 

times<-1:45

dat<-as.matrix(read.csv("neon_cx_tarsalis.csv", header = T))
dat2<-dat[,-1] ## removes tBin
dat2<-dat2[,-3] ## removes JORN
dat2<-t(dat2)

x<-cleandat(dat2,times,3)$cdat

## Get synchrony preserving aaft surrogates of your mosquito data
nsurrogs <- 1000 # for a quick test 
surrtype <- "aaft"
syncpres <- TRUE
sur <- surrog(x,nsurrogs,surrtype,syncpres)

df_slope <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(df_slope) <- c("slope_surrog","slope_surrog2")


# Define a function to process each iteration
process_iteration <- function(i) {
  x <- i  # Surrogate data
  
  moswpmf <- wpmf(x, times = times, sigmethod = "aaft")
  
  g <- abs(moswpmf$value)
  g12 <- g[12, -46]
  g8_14 <- colMeans(g[8:14, -46], na.rm = TRUE)
  
  out1 <- lm(g12 ~ times)
  slope_surrog <- out1$coefficients[2]
  names(slope_surrog) <- "slope_surrog"
  
  out2 <- lm(g8_14 ~ times)
  slope_surrog2 <- out2$coefficients[2]
  names(slope_surrog2) <- "slope_surrog2"
  
  # Return the calculated slopes for this iteration
  return(c(slope_surrog, slope_surrog2))
}



############# try optimize run time for getting surrogate data (to be run on hpg) #############

# ## Parallel processing in principle:
# # start parallel processing using doParallel package
# 
# ncores <- detectCores(logical = FALSE) # does not include 'hyperthreaded cores'
# cl <- makeCluster(floor(0.25 * ncores))
# registerDoParallel(cl)
# 
# # stop cluster explicitly using stopCluster to free up the computational resources
# stopCluster(cl) # clusters automatically stop at the end of the R session



#Set the number of cores for parallel processing
num_cores <- 20  # Adjust as needed
#num_cores <- detectCores(logical = FALSE)

# Create a parallel cluster and register it
cl <- makeCluster(num_cores)
#cl <- makeCluster(floor(0.75 * num_cores)
registerDoParallel(cl)

# Apply the function in parallel
results <- foreach(j = 1:length(sur), .combine = rbind, .packages = "wsyn") %dopar% {
  process_iteration(sur[[j]])
}

# Stop the parallel cluster
stopCluster(cl)



# Combine the results into a data frame
df_slope <- as.data.frame(results)

# Save the resulting slope data frame
file.name <- str_c("Surrogate_slopes_", surrtype, "_", length(sur), "_", hms, ".csv")
write.csv(df_slope, file = file.name)



