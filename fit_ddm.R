fit_ddm <- function(data, rt, response, startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)) {
    
    # install packages
    packages <- c("rtdists", "ucminf")
    toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(toInstall)) install.packages(toInstall)
    rm(packages); rm(toInstall) 
    library(rtdists); library(ucminf)
    
    # if response coded as 0 or 1, recode as 'lower' and 'upper'
    data[data[, response] == 0, response] <- 'lower'
    data[data[, response] == 1, response] <- 'upper'
    
    # calculate likelihood of model parameters given observed data
    likelihood_ddm <- function(params, rt, response) {
        densities <- ddiffusion(rt = rt, response = response, a = params['a'], v = params['v'], t0 = params['t0'], z = params['z'] * params['a']) + .Machine$double.eps
        return(-2*sum(log(densities))) # return -2 * log-likelihood
    }
    
    # startParams <- c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)
    # nlminb(startParams, likelihood_ddm, rt = data$rt, response = data$response)
    optimResults <- ucminf(startParams, likelihood_ddm, rt = data[, rt], response = data[, response])
    return(optimResults)
}

# #### test function
# source('https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ddm.R') # function also on github
# data <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2)
# results_ddm <- fit_ddm(data = data, rt = 'rt', response = 'response') # response argument accepts values coded as 0/1 or 'lower'/'upper'
# results_ddm$par
# 
# # compare with fit_ezddm
# source("https://raw.githubusercontent.com/hauselin/Rcode/master/fit_ezddm.R")
# data$response2 <- ifelse(data$response == 'upper', 1, 0) # recode response to 0 and 1
# results_ezddm <- fit_ezddm(data = data, reactiontime = 'rt', acc = 'response2')
# results_ezddm[, 1:2] <- results_ezddm[, 1:2] * 10 # multiply by ten to match fit_ddm results
# results_ezddm