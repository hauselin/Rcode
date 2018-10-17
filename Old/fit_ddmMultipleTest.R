library(data.table); library(rtdists)
data1 <- rdiffusion(n = 300, a = 2, v = 0.3, t0 = 0.5, z = 0.3 * 2)
data2 <- rdiffusion(n = 300, a = 2, v = -0.3, t0 = 0.5, z = 0.3 * 2)
setDT(data1)
setDT(data2)
data1[, id := 1]
data2[, id := 2]
data <- bind_rows(data1, data2)
data[, cond1 := sample(c("a", "b"), 600, replace = T)]
data[, cond2 := sample(c("c", "d"), 600, replace = T)]

startParams = c(a = 2, v = 0.1, t0 = 0.3, z = 0.5)

# define function to calculate likelihood of model parameters given observed data
likelihood_ddm <- function(params, rt, response) {
    densities <- ddiffusion(rt = rt, response = response, 
                            a = params['a'], 
                            v = params['v'], 
                            t0 = params['t0'], 
                            z = params['z'] * params['a']) + .Machine$double.eps
    return(-2*sum(log(densities))) # return -2 * sum(log-likelihood)
}



nlminb(startParams, likelihood_ddm, rt = data$rt, response = data$response)
ucminf(startParams, likelihood_ddm, rt = data$rt, response = data$response)

data[, nlminb(startParams, likelihood_ddm, rt = rt, response = response)]
res <- data[, nlminb(startParams, likelihood_ddm, rt = rt, response = response), by = .(id, cond1, cond2)]
res2 <- data[, ucminf(startParams, likelihood_ddm, rt = rt, response = response), by = .(id, cond1, cond2)]
res[, parName := c("a", "v", "t0", "z")]

data[, nlminb(startParams, likelihood_ddm, rt = rt, response = response), by = .(id)]
data[, ucminf(startParams, likelihood_ddm, rt = rt, response = response), by = .(id)]

data[, ucminf(startParams, likelihood_ddm, rt = rt, response = response)[c('par', 'value', 'convergence')], by = .(id)]


formulaString <- paste0(id)
if (!is.null(group)) {
    for (i in 1:length(group)) {
        formulaString <- paste(formulaString, group[i], sep = " + ")
    }
    formulaString <- paste(formulaString, "iterations", "convergence", "objective", sep = " + ")
} 
formulaString <- paste0(formulaString, " ~ parName")
form <- as.formula(formulaString)

resultsWide <- dcast(data = res, formula = form, value.var = c('par'))

resultsWide2 <- left_join(trials, resultsWide)
resultsFinal <- select(resultsWide2, id, group, N, a:z, everything())

# remove temporary_subject variable
if (id == 'temporary_subject') {
    resultsFinal$temporary_subject <- NULL
}



fit_ddm(data = data, rts = "rt", responses = "response", id = "id")
fit_ddm(data = data, rts = "rt", responses = "response", id = "id", group = "cond1")
fit_ddm(data = data, rts = "rt", responses = "response", id = "id", group = c("cond1", "cond2"))
fit_ddm(data = data, rts = "rt", responses = "response")

