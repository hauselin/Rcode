# install really useful packages
packages <- c("dplyr", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) {
    install.packages(toInstall)
} else {
    library(dplyr); library(data.table)
}
rm(packages); rm(toInstall)

fit_ezddm <- function(data, reactiontime, accuracy, id = NULL, group = NULL) {
    
    message("Reaction times must be in seconds.\nAccuracy or choice must be coded as 0 or 1.")
    
    setDT(data) # convert to data table
    
    # if no id variable provided, assume it's just one subject's data
    if (is.null(id)) {
        id <- "temporary_subject"
        data[, (id) := 1] 
        message("id variable not provided. Assuming single-subject data.")
    }
    
    # for accurate responses (coded as 1), calculate mean RT and RT variance for each subject, each condition
    ddmRt <- data[get(accuracy) == 1, .(rt = mean(get(reactiontime), na.rm = T), rtVar = var(get(reactiontime), na.rm = T)), by = c(id, group)]
    
    # calculate accuracy for each subject, each condition
    ddmAcc <- data[, .(acc = mean(get(accuracy), na.rm = T), n = .N), by = c(id, group)]
    
    if (sum(ddmAcc[, acc] %in% c(0, 0.5, 1)) > 0) {
        n_corrected <- sum(ddmAcc[, acc] %in% c(0, 0.5, 1))
        message(paste0("Mean accuracies (n = ", n_corrected, ") that are 0, 0.5, or 1 have been adjusted slightly for model fitting."))
        ddmAcc[, acc_adjust := 0]
        ddmAcc[acc %in% c(0, 0.5, 1), acc_adjust := 1]
    }
    
    # if acc is 1, apply edge correction
    ddmAcc[acc == 1, acc := edgeCorrect(n)] # edge correction
    # if acc is 0 or 50, add 0.001 to acc a bit so model fitting works
    ddmAcc[acc %in% c(0, 0.5), acc := acc + 0.0001]
    
    dataForDDM <- left_join(ddmRt, ddmAcc, by = c(id, group))
    
    # fit ez ddm model to each subject, each condition
    ddmResults <- dataForDDM[, ezddm(propCorrect = acc, rtCorrectVariance_seconds = rtVar, rtCorrectMean_seconds = rt), by = c(id, group)]
    
    ddmResults <- left_join(ddmResults, ddmRt, by = c(id, group)) %>% left_join(ddmAcc, by = c(id, group))
    
    # remove temporary_subject variable
    if (id == 'temporary_subject') {
        ddmResults[, temporary_subject := NULL]
    }
    
    setDT(ddmResults) # ensure it's data table format
    setnames(ddmResults, c("v", "a", "Ter", "rt", "rtVar"), c("drift_v", "threshold_a", "ndt_Ter", "rt_acc1", "rtVar_acc1"))
    return(ddmResults[])
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
            cat("Oops, propCorrect == 0\n")
        } else if (propCorrect == 0.5) {
            cat("Oops, propCorrect == 0.5\n")
        } else if (propCorrect == 1) {
            cat("Oops, propCorrect == 1\n")
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
