# install really useful packages
packages <- c("dplyr", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) {
    install.packages(toInstall)
} else {
    library(dplyr); library(data.table); library(piecewiseSEM); library(compute.es); library(sjstats)
}
rm(packages); rm(toInstall)

fit_ezddm <- function(data, rtvar, accvar, idvar, groupvar) {
    
    ddmRt <- data[get(accvar) == 1, .(rt = mean(get(rtvar), na.rm = T), rtVar = var(get(rtvar), na.rm = T)), 
                  by = mget(c(idvar, groupvar))]
    
    ddmAcc <- data[, .(acc = mean(get(accvar), na.rm = T), n = .N), by = mget(c(idvar, groupvar))]
    ddmAcc[acc == 1, acc := edgeCorrect(n)]
    
    dataForDDM <- left_join(ddmRt, ddmAcc)
    
    ddmResults <- dataForDDM[, ezddm(propCorrect = acc, rtVar, rtCorrectMean_seconds = rt), by = mget(c(idvar, groupvar))]
    
    ddmResults <- left_join(ddmResults, ddmRt) %>% left_join(ddmAcc)
    
    setDT(ddmResults)
    return(ddmResults)
}

ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds) {
    
    #'  propCorrect: proportion correct (apply edge correction if necessary)
    #'  rtVariance: variance of correct reaction times (in seconds)
    #'  rtMean: mean of correct reaction times (in seconds)
    
    s <- 0.1 # s is scaling parameter (defaults to 0.1 in Ratcliff's models)
    s2 <- s^2 # variance
    
    # if propCorrect equals 0, 0.5, or 1, this method will not work, and an edge correction is required
    if (propCorrect %in% c(0, 0.5, 1)) {
        
        v <- as.numeric(NA)
        a <- as.numeric(NA)
        Ter <- as.numeric(NA)
        
        if (propCorrect == 0) {
            cat("Oops, propCorrect == 0!\n")
        } else if (propCorrect == 0.5) {
            cat("Oops, propCorrect == 0.5!\n")
        } else if (propCorrect == 1) {
            cat("Oops, propCorrect == 1!\n")
        }
        
    } else {
        
        L <- qlogis(propCorrect) # calculates logit
        x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
        v <- sign(propCorrect - 0.5) * s * x^(1/4) # drift rate
        a <- s2 * qlogis(propCorrect)/v # threshold 
        y <- -v*a/s2
        MDT <- (a/(2*v)) * (1-exp(y))/(1 + exp(y))
        Ter <- rtCorrectMean_seconds - MDT # non-decision time
        
    }
    
    return(data.frame(v, a, Ter))
}

# test function
# ezddm(.802, .112, .723)
# ezddm(0.8881988, 0.1005484, 0.9010186)
# library(EZ2)
# Data2EZ(.802, .112, .723)
# Data2EZ(0.8881988, 0.1005484, 0.9010186)
# data.frame(Data2EZ(.802, .112, .723))

edgeCorrect <- function(n) {
    #' n: number of observations
    return(1 - (1 / (2 * n)))
}
