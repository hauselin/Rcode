# install packages
packages <- c("rtdists", "ucminf", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
rm(packages); rm(toInstall) 
library(rtdists); library(ucminf); library(data.table)

fit_ddm <- function(data, rts, responses, startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)) {
    
    setDT(data)
    
    # check if rt is in seconds
    if (data[, mean(get(rts), na.rm = T)] > 10) {
        message("Check if reaction time is in seconds, not milliseconds!")
    } 
    
    # remove rts or responses rows
    data <- data[!is.na(get(rts)), ]
    data <- data[!is.na(get(responses)), ]
    
    # if response coded as 0 or 1, recode as 'lower' and 'upper'
    if (data[, unique(get(responses))][1] %in% c(0, 1)) {
        data[, response := as.character(get(responses))]
        data[, response := ifelse(get(responses) == "1", 'upper', 'lower')]
    }

    # calculate likelihood of model parameters given observed data
    likelihood_ddm <- function(params, rt, response) {
        densities <- ddiffusion(rt = rt, response = response, 
                                a = params['a'], 
                                v = params['v'], 
                                t0 = params['t0'], 
                                z = params['z'] * params['a']) + .Machine$double.eps
        return(-2*sum(log(densities))) # return -2 * sum(log-likelihood)
    }
    
    # startParams <- c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)
    # optimResults <- nlminb(startParams, likelihood_ddm, rt = data[, get(rts)], response = data[, get(responses)])
    optimResults <- ucminf(startParams, likelihood_ddm, rt = data[, get(rts)], response = data[, get(responses)])
    return(optimResults)
}

#### test function
source('https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ddm.R') # function also on github
data <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2)
results_ddm <- fit_ddm(data = data, rts = 'rt', responses = 'response') # response argument accepts values coded as 0/1 or 'lower'/'upper'
results_ddm$par

data <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2)
rts <- "rt"
responses <- "response"
data$response <- ifelse(data$response == 'upper', 1, 0) # recode response to 0 and 1
setDT(data)
data[c(1, 3), rt := NA]
startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)


check <- rdiffusion(n = 300, a = results_ddm$par['a'], v = results_ddm$par['v'], t0 = results_ddm$par['t0'], z = results_ddm$par['z'])
tbl_dt(data)[, mean(rt), response]
tbl_dt(check)[, mean(rt), response]

# compare with fit_ezddm
source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ezddm.R")
data$response2 <- ifelse(data$response == 'upper', 1, 0) # recode response to 0 and 1
results_ezddm <- fit_ezddm(data = data, reactiontime = 'rt', acc = 'response2')
results_ezddm[, 1:2] <- results_ezddm[, 1:2] * 10 # multiply by ten to match fit_ddm results
results_ezddm