# install really useful packages
packages <- c("dplyr", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) {
    install.packages(toInstall)
} else {
    library(dplyr); library(data.table)
}
rm(packages); rm(toInstall)

fit_ezddm <- function(data, rt, acc, id, group) {
    
    message("Reaction times (rt) must be in seconds.\n
        Accuracy or choice (acc) must be coded as 0 or 1.")
    
    setDT(data) # convert to data table
    
    # for accurate responses (coded as 1), calculate mean RT and RT variance for each subject, each condition
    ddmRt <- data[get(acc) == 1, .(rt = mean(get(rt), na.rm = T), rtVar = var(get(rt), na.rm = T)), by = c(idvar, group)]
    
    # calculate accuracy for each subject, each condition
    ddmAcc <- data[, .(acc = mean(get(acc), na.rm = T), n = .N), by = c(id, group)]
    ddmAcc[acc == 1, acc := edgeCorrect(n)] # edge correction
    
    dataForDDM <- left_join(ddmRt, ddmAcc)
    
    # fit ez ddm model to each subject, each condition
    ddmResults <- dataForDDM[, ezddm(propCorrect = acc, rtVar, rtCorrectMean_seconds = rt), by = c(id, group)]
    
    ddmResults <- left_join(ddmResults, ddmRt) %>% left_join(ddmAcc)
    
    setDT(ddmResults) # ensure it's data table format
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
