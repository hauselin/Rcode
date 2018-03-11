reportCHISQ <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    pdigits <- paste0("%.", decimal + 1, "f") # p values always 1 more digit
    
    # example output: r(30) = 0.82, p < .001
    estimates <- data.frame(df = model$parameter, 
                            statistic = model$statistic,
                            p.value = model$p.value)
    rownames(estimates) <- NULL
    
    # effect sizes
    computeES <- chies(chi.sq = estimates$statistic, n = sum(model$observed), verbose = F)
    estimates$es.r <- computeES$r
    
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
    
    formattedOutput <- paste0("X2(", estimatesRound$df, ")",
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