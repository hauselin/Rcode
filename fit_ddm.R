# install packages
packages <- c("rtdists", "ucminf", "data.table", "tidyverse", "dtplyr", "doFuture")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
rm(packages); rm(toInstall) 
library(rtdists); library(ucminf); library(data.table); library(tidyverse); library(dtplyr); library(doFuture)

fit_ddm <- function(data, rts, responses, id = NULL, group = NULL, startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5), simCheck = TRUE, decimal = 4, parallel = FALSE) {
    
    if (parallel) {
        registerDoFuture()
        plan(multiprocess)
    }
    
    data <- tbl_dt(data)
    
    # create new variables
    data$rtCol <- data[, get(rts)]
    data$responseCol <- data[, get(responses)]
    
    # check if rt is in seconds
    if (mean(data$rtCol, na.rm = T) > 100) {
        stop("Check if reaction time is in seconds, not milliseconds!")
    }
    
    # if no id variable provided, assume it's just one subject's data
    if (is.null(id)) {
        id <- "temporary_subject"
        data$temporary_subject <- 1
    }
    
    # remove rts or responses rows
    data <- data[!is.na(rtCol), ]
    data <- data[!is.na(responseCol), ]
    
    # recode response accordingly (character and integer)
    if (data[, unique(responseCol)][1] %in% c(0, 1)) {
        data[, response_num := responseCol]
        data[, response_char := ifelse(response_num == 1, 'upper', 'lower')]
    } else if (data[, unique(responseCol)][1] %in% c('upper', 'lower')) {
        data[, response_char := responseCol]
        data[, response_num := ifelse(response_char == "upper", 1, 0)]
    }
    
    # compute rt and responses
    behavOverall <- data[, .(response = round(mean(response_num, na.rm = T), 3), rtOverall = round(mean(rtCol, na.rm = T), 3)), by = c(id, group)]
    behav0 <- data[response_num == 0, .(rt0 = round(mean(rtCol, na.rm = T), 3)), by = c(id, group)]
    behav1 <- data[response_num == 1, .(rt1 = round(mean(rtCol, na.rm = T), 3)), by = c(id, group)]
    behav <- left_join(behavOverall, behav0, by = c(id, group)) %>% left_join(behav1, by = c(id, group))
    
    # define function to calculate likelihood of model parameters given observed data
    likelihood_ddm <- function(params, rt, response) {
        if (params['t0'] < 0.03 | params['a'] <= 0 | params['z'] <= 0 | params['z'] > 1) return(1e6)
        densities <- ddiffusion(rt = rt, response = response, 
                                a = params['a'], 
                                v = params['v'], 
                                t0 = params['t0'], 
                                z = params['z'] * params['a']) + .Machine$double.eps
        return(-2*sum(log(densities))) # return -2 * sum(log-likelihood)
    }
    
    dataGroup <- data[, .(n = .N), by = c(id, group)]
    dataGroup0 <- data[response_num == 0, .(n0 = .N), by = c(id, group)]
    dataGroup1 <- data[response_num == 1, .(n1 = .N), by = c(id, group)]
    dataGroup <- left_join(dataGroup, dataGroup0, by = c(id, group)) %>% left_join(dataGroup1, by = c(id, group))
    
    # startParams <- c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)
    # optimResults <- nlminb(startParams, likelihood_ddm, rt = data[, get(rts)], response = data[, get(responses)])
    # optimResults <- ucminf(startParams, likelihood_ddm, rt = data[, get(rts)], response = data[, get(responses)])
    # return(optimResults)
    # optimize for each subject, each condition/group
    # res <- data[, nlminb(startParams, likelihood_ddm, rt = get(rts), response = get(responses)), by = c(id, group)] # nlminb optimization
    if (any(class(startParams) %in% c("data.frame"))) {
        
        if (parallel) {
            message("Running parallel loops for each combination of starting values...")
            res <- foreach(startParamsI = 1:nrow(startParams)) %dopar% {
                # startParametersTemp <- c(a = startParams$a[startParamsI], v = startParams$v[startParamsI], t0 = startParams$t0[startParamsI], z = startParams$z[startParamsI])
                data[, ucminf(c(a = startParams$a[startParamsI], v = startParams$v[startParamsI], t0 = startParams$t0[startParamsI], z = startParams$z[startParamsI]), likelihood_ddm, rt = rtCol, response = response_char)[c('par', 'value', 'convergence')], by = c(id, group)] # ucminf optimization
            }
            res <- rbindlist(res)
            
        } else {
            message("Starting parameters for optimization:")
            res <- data.frame()
            for (startParamsI in 1:nrow(startParams)) {
                startParametersTemp <- c(a = startParams$a[startParamsI], v = startParams$v[startParamsI], t0 = startParams$t0[startParamsI], z = startParams$z[startParamsI])
                print(startParametersTemp)
                resTemp <- data[, ucminf(startParametersTemp, likelihood_ddm, rt = rtCol, response = response_char)[c('par', 'value', 'convergence')], by = c(id, group)] # ucminf optimization
                res <- bind_rows(res, resTemp)
            }
        }
        
        setDT(res)
        res <- distinct(res) # get rid of results/rows with exactly same results
        colNamesOriginal <- names(res) # save column order
        res <- left_join(res[, .(value = min(value)), by = c(id, group)], res, by = c(id, group, "value")) # find minimum values by group
        setcolorder(res, names(colNamesOriginal)) # reorder columns to original order
            
    } else {
        message("Starting parameters for optimization:")
        print(startParams)
        res <- data[, ucminf(startParams, likelihood_ddm, rt = rtCol, response = response_char)[c('par', 'value', 'convergence')], by = c(id, group)] # ucminf optimization
        setDT(res)
    }
    res[, parName := c("a", "v", "t0", "z")]
    
    # convert long to wide form
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
    resultsFinal <- select(resultsWide2, id, group, n:n1, a, v, t0, z, everything())
    resultsFinal <- left_join(resultsFinal, behav, by = c(id, group))
    setDT(resultsFinal)
    # print(resultsFinal)
    
    # simulate data
    if (simCheck) {
        resultsToSimulate <- copy(resultsFinal)
        simulatedData <- resultsToSimulate[, rdiffusion(n = 1000, a = a, v = v, t0 = t0, z = z * a), by = c(id, group)]
        simulatedData[, response_num := as.numeric()]
        simulatedData[response == 'upper', response_num := 1]
        simulatedData[response == 'lower', response_num := 0]
        simulateBehavOverall <- simulatedData[, .(responseSim = round(mean(response_num, na.rm = T), 3), rtOverallSim = round(mean(rt, na.rm = T), 3)), by = c(id, group)]
        simulateBehav0 <- simulatedData[response_num == 0, .(rt0Sim = round(mean(rt, na.rm = T), 3)), by = c(id, group)]
        simulateBehav1 <- simulatedData[response_num == 1, .(rt1Sim = round(mean(rt, na.rm = T), 3)), by = c(id, group)]
        simulateBehav <- left_join(simulateBehavOverall, simulateBehav0, by = c(id, group)) %>% left_join(simulateBehav1, by = c(id, group))
        resultsFinal <- left_join(resultsFinal, simulateBehav, by = c(id, group))
        resultsFinal <- select(resultsFinal, 1, group, n:value, response, responseSim, rtOverall, rtOverallSim, rt0, rt0Sim, rt1, rt1Sim)
    }
    
    # round results
    resultsFinal[, a := round(a, decimal)]
    resultsFinal[, v := round(v, decimal)]
    resultsFinal[, t0 := round(t0, decimal)]
    resultsFinal[, z := round(z, decimal)]
    
    # remove temporary_subject variable
    if (id == "temporary_subject") {
        resultsFinal[, temporary_subject := NULL]
    }
    
    return(tbl_dt(resultsFinal))
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