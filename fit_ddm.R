# install packages
packages <- c("rtdists", "ucminf", "data.table", "tidyverse")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
rm(packages); rm(toInstall) 
library(rtdists); library(ucminf); library(data.table); library(tidyverse)

fit_ddm <- function(data, rts, responses, id = NULL, group = NULL, startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)) {
    
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
    
    # if no id variable provided, assume it's just one subject's data
    if (is.null(id)) {
        id <- "temporary_subject"
        data[, (id) := 1] 
        # message("id variable not provided. Assuming single-subject data.")
    }

    # define function to calculate likelihood of model parameters given observed data
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
    # optimResults <- ucminf(startParams, likelihood_ddm, rt = data[, get(rts)], response = data[, get(responses)])
    # return(optimResults)
    
    dataGroup <- data[, .(n = .N), by = c(id, group)]
    
    # optimize for each subject, each condition/group
    # res <- data[, nlminb(startParams, likelihood_ddm, rt = get(rts), response = get(responses)), by = c(id, group)] # nlminb optimization
    res <- data[, ucminf(startParams, likelihood_ddm, rt = get(rts), response = get(responses))[c('par', 'value', 'convergence')], by = c(id, group)] # ucminf optimization
    res[, parName := c("a", "v", "t0", "z")]
    
    formulaString <- paste0(id)
    if (!is.null(group)) {
        for (i in 1:length(group)) {
            formulaString <- paste(formulaString, group[i], sep = " + ")
        }
    } 
    formulaString <- paste(formulaString, "convergence", "value", sep = " + ")
    formulaString <- paste0(formulaString, " ~ parName")
    form <- as.formula(formulaString)
    # print(form)
    
    resultsWide <- dcast(data = res, formula = form, value.var = c('par'))
    
    resultsWide2 <- left_join(dataGroup, resultsWide, by = c(id, group))
    resultsFinal <- select(resultsWide2, id, group, n, a, v, t0, z, everything())
    # print(resultsFinal)
    
    # remove temporary_subject variable
    if (id == 'temporary_subject') {
        resultsFinal$temporary_subject <- NULL
    }
    
    setDT(resultsFinal)
    return(resultsFinal[])
}

#### test function
# library(data.table); library(rtdists); 
# data1 <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2) # simulate data
# data2 <- rdiffusion(n = 300, a = 2, v = -0.3, t0 = 0.5, z = 0.3 * 2) # simulate data
# setDT(data1) # convert to data.table
# setDT(data2)
# data1[, subject := 1] # add subject id
# data2[, subject := 2]
# dataAll <- rbind(data1, data2)
# dataAll[, cond1 := sample(c("a", "b"), 600, replace = T)] # randomly assign conditions a/b
# dataAll[, cond2 := sample(c("c", "d"), 600, replace = T)] # randomly assign conditions c/d
# fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject") # fit model to each subject (no conditions)
# fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1") # fit model to each subject by cond1
# fit_ddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2")) # fit model to each subject by cond1,cond2
# fit_ddm(data = dataAll, rts = "rt", responses = "response") # fit model to just entire data set (1 subject, 1 condition) 