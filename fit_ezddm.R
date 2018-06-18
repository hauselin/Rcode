# fit_ezddm function fits Wagenmaker et al.'s (2007) EZ-diffusion model for two-choice response time tasks. To use the function, ensure your dataframe is in long form, has single-trial reaction time (in seconds) and responses (coded as 0 or 1) on each row. You can use the function to fit the EZ-diffusion model to just a single subject or multiple subjects, and separately for each experimental condition (see below for examples).

# install really useful packages
packages <- c("dplyr", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
rm(packages); rm(toInstall)
library(dplyr); library(data.table)

fit_ezddm <- function(data, rts, responses, id = NULL, group = NULL) {
    
    # message("Fits EZ-diffusion model (Wagenmaker et al., 2007, Psychonomic Bulletin & Review).\nResponses or choice must be coded as 0 (lower bound) or 1 (upper bound).")
    
    setDT(data) # convert to data table
    
    if (data[, mean(get(rts), na.rm = T)] > 10) {
        message("Check if reaction time is in seconds, not milliseconds!")
    }
    
    if (class(data[1, get(responses)]) %in% c('character', 'factor')) {
        stop("Check if responses are coded as 0 (lower) and 1 (upper)!")
    }
    
    # if no id variable provided, assume it's just one subject's data
    if (is.null(id)) {
        id <- "temporary_subject"
        data[, (id) := 1] 
    }
    
    # get grouping variables
    dataGroup <- data[, .(trials = .N), by = c(id, group)]
    dataGroup[, trials := NULL]
    
    # for accurate responses (coded as 1), calculate mean RT and RT variance for each subject, each condition
    ddmRt <- data[get(responses) == 1, 
                  .(rt = mean(get(rts), na.rm = T), rtVar = var(get(rts), na.rm = T)),
                  by = c(id, group)]
    
    # calculate responses for each subject, each condition
    ddmAcc <- data[, .(acc = mean(get(responses), na.rm = T), n = .N), by = c(id, group)]
    
    if (sum(ddmAcc[, acc] %in% c(0.5, 1)) > 0) {
        n_corrected <- sum(ddmAcc[, acc] %in% c(0.5, 1))
        message(paste0("Mean accuracies (n = ", n_corrected, ") that are 0.5, or 1 have been adjusted slightly for model fitting."))
        ddmAcc[, acc_adjust := 0]
        ddmAcc[acc %in% c(0.5, 1), acc_adjust := 1]
    }
    
    # if acc is 1, apply edge correction
    ddmAcc[acc == 1, acc := edgeCorrect(n)] # edge correction
    # if acc is 0 or 50, add 0.001 to acc a bit so model fitting works
    ddmAcc[acc %in% c(0.5), acc := acc + 0.00001]
    
    dataForDDM <- left_join(ddmRt, ddmAcc, by = c(id, group))
    setDT(dataForDDM)
    
    # fit ez ddm model to each subject, each condition
    ddmResults <- dataForDDM[, ezddm(propCorrect = acc, rtCorrectVariance_seconds = rtVar, rtCorrectMean_seconds = rt), by = c(id, group)]
    
    ddmResults <- left_join(dataGroup, ddmResults, by = c(id, group)) %>% 
        left_join(ddmRt, by = c(id, group)) %>%
        left_join(ddmAcc, by = c(id, group))
    
    ddmResults <- select(ddmResults, id, group, n, everything()) # reorder columns
    # remove temporary_subject variable
    if (id == 'temporary_subject') {
        ddmResults$temporary_subject <- NULL
    }
    setDT(ddmResults) # ensure it's data table format
    setnames(ddmResults, c("a", "v", "Ter", "rt", "rtVar"), c("a_threshold", "v_drift", "ndt_Ter", "rt_correct", "rtVar_correct"))
    return(ddmResults[])
}

ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds, nTrials = NULL) {
    
    #'  propCorrect: proportion correct (apply edge correction if necessary)
    #'  rtVariance: variance of correct reaction times (in seconds)
    #'  rtMean: mean of correct reaction times (in seconds)
    #'  nTrials (optional): number of trials (useful for edge correction)
    
    s <- 0.1 # s is scaling parameter (defaults to 0.1 in Ratcliff's models)
    s2 <- s^2 # variance
    
    v <- as.numeric(NA)
    a <- as.numeric(NA)
    Ter <- as.numeric(NA)
    
    # if propCorrect equals 0, 0.5, or 1, this method will not work, and an edge correction is required
    if (propCorrect %in% c(0, 0.5, 1)) {
        
        if (propCorrect == 0) {
            return(cat("Oops, propCorrect == 0. Can't fit model! D:"))
        } else if (propCorrect == 0.5) {
            cat("Oops, propCorrect == 0.5 (chance performance; drift will be close to 0). Added 0.00001 to propCorrect.\n")
            propCorrect <- propCorrect + 0.00001
        } else if (propCorrect == 1) {
            if (!is.null(nTrials)) {
                cat("Oops, propCorrect == 1. Applied edge correction.\n")
                propCorrect <- 1 - (1 / (2 * nTrials))
            } else {
                cat("Oops, propCorrect == 1. Edge correction required. Provide number of trials (nTrials).\n")    
            }
        }
        
    } 
    
    if (propCorrect != 1) {
    
        L <- qlogis(propCorrect) # calculates logit
        x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
        v <- sign(propCorrect - 0.5) * s * x^(1/4) # drift rate
        a <- s2 * qlogis(propCorrect)/v # threshold 
        y <- -v*a/s2
        MDT <- (a/(2*v)) * (1-exp(y))/(1 + exp(y))
        Ter <- rtCorrectMean_seconds - MDT # non-decision time

    }
    
    return(round(data.frame(a, v, Ter), 6))
}

edgeCorrect <- function(n) {
    #' n: number of observations
    return(1 - (1 / (2 * n)))
}


# test function
# library(rtdists)
# rt1 <- rdiffusion(200, a=3, v=2, t0=0.5)
# rt1$response <- ifelse(rt1$response == "upper", 1, 0)
# rt1$rt2 <- rt1$rt * 1000
# fit_ezddm(data = rt1, rts = "rt", responses = "response")
# fit_ezddm(data = rt1, rts = "rt2", responses = "response")


data1 <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2) # simulate data
data2 <- rdiffusion(n = 300, a = 2, v = -0.3, t0 = 0.5, z = 0.3 * 2) # simulate data
setDT(data1) # convert to data.table
setDT(data2)
data1[, subject := 1] # add subject id
data2[, subject := 2]
dataAll <- rbind(data1, data2)
dataAll[, cond1 := sample(c("a", "b"), 600, replace = T)] # randomly assign conditions a/b
dataAll[, cond2 := sample(c("c", "d"), 600, replace = T)] # randomly assign conditions c/d
dataAll$response <- ifelse(dataAll$response == "upper", 1, 0)
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject") # fit model to each subject (no conditions)
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1") # fit model to each subject by cond1
fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2")) # fit model to each subject by cond1,cond2
fit_ezddm(data = dataAll, rts = "rt", responses = "response") # fit model to just entire data set (1 subject, 1 condition) 

# ezddm(.802, .112, .723)
# ezddm(.5, .112, .723)
# ezddm(.51, .112, .723)
# ezddm(0, .112, .723)
# ezddm(0.0001, .112, .723)
# ezddm(0, .112, .723)
# ezddm(0.005, .112, .723)
# ezddm(0.005, .112, .723)
# ezddm(1, .112, .723, 100)
#
# ezddm(0.8881988, 0.1005484, 0.9010186)
# library(EZ2)
# Data2EZ(.802, .112, .723)
# Data2EZ(.5, .112, .723)
# Data2EZ(0.8881988, 0.1005484, 0.9010186)
# Data2EZ(0.1, 0.1005484, 0.9010186)
# Data2EZ(0.00001, 0.1005484, 0.9010186)
# ezddm(0.000001, 0.1005484, 0.9010186)
# ezddm(0.00001, 0.1005484, 0.9010186)
# ezddm(0.5, 0.1005484, 0.9010186)
# ezddm(0.51, 0.1005484, 0.9010186)
# data.frame(Data2EZ(.802, .112, .723))