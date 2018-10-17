reportTtest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    pdigits <- paste0("%.", decimal + 1, "f") # p values always 1 more digit
    
    # example output: t(30) = 5.82, p < .001
    estimates <- data.frame(df = model$parameter, 
                            statistic = model$statistic,
                            p.value = model$p.value)
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # make a copy of estimates and convert to correct dp
    estimatesRound <- sapply(estimates, function(x) ifelse(abs(x) < 0.01, round(x, 3), round(x, decimal)))
    estimatesRound <- sapply(estimatesRound, function(x) ifelse(abs(x) < 0.01, sprintf(pdigits, x), sprintf(digits, x)))
    estimatesRound <- data.frame(as.list(estimatesRound), stringsAsFactors = FALSE)
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    
    formattedOutput <- paste0("t(", estimatesRound$df, ")",
                              " = ", estimatesRound$statistic, 
                              ", p ", estimatesRound$p.value, 
                              ", r = ", estimatesRound$es.r)
    
    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "−", formattedOutput)
    }
    
    formattedOutputDf <- data.frame(results = as.character(formattedOutput), stringsAsFactors = FALSE)
    
    outputList <- list(results = formattedOutputDf)
    
    if (showTable) {
        
        # format table nicely
        estimatesOutput <- apply(estimates, 2, round, decimal + 1)
        estimatesOutput <- data.frame(term = as.character(model$data.name), as.list(estimatesOutput), stringsAsFactors = FALSE)
        outputList$results2 <- estimatesOutput
    }
    
    if (showEffectSizesTable) {
        
        # get all other effect sizes
        effectSizes <- es(r = abs(round(as.numeric(estimates$es.r), decimal + 1)), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(model$data.name), as.list(effectSizes), stringsAsFactors = FALSE)
        
    }
    
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}

#### examples ####
# model <- t.test(mtcars$mpg, mtcars$vs)
# reportTtest(model)
# reportTtest(model, showTable = T, decimal = 3)
# reportTtest(model, showTable = T, showEffectSizesTable = T, decimal = 4)
# reportTtest(model, showTable = F, showEffectSizesTable = T, decimal = 2)