reportLM <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    pdigits <- paste0("%.", decimal + 1, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    
    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- df.residual(model) # get model degrees of freedom
    estimates <- estimates[, c(1, 2, 5, 3, 4)] # sort columns
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # R2
    estimates$es.r.squared <- summary(model)$r.squared
    estimates$es.adj.r.squared <- summary(model)$adj.r.squared
    
    # make a copy of estimates and convert o correct dp
    estimatesRound <- apply(estimates[, -1], 2, function(x) ifelse(abs(x) < 0.01, round(x, 3), round(x, decimal)))
    estimatesRound <- apply(estimatesRound, 2, function(x) ifelse(abs(x) < 0.01, sprintf(pdigits, x), sprintf(digits, x)))
    estimatesRound <- data.frame(term = estimates$term, estimatesRound, stringsAsFactors = FALSE)
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    
    formattedOutput <- paste0("b = ", estimatesRound$estimate, 
                              ", SE = ", estimatesRound$std.error, 
                              ", t(", estimatesRound$df, ")", 
                              " = ", estimatesRound$statistic,
                              ", p ", estimatesRound$p.value, 
                              ", r = ", estimatesRound$es.r)
    
    # convert hyphens to minus (only possible on UNIX systems)
    if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative sign is minus, not hyphens
        formattedOutput <- gsub("-", replacement = "−", formattedOutput)
    }
    
    formattedOutputDf <- data.frame(term = as.character(estimates$term), 
                                    results = as.character(formattedOutput),
                                    stringsAsFactors = FALSE)
    
    outputList <- list(results = formattedOutputDf)
    
    if (showTable) {
        
        # format table nicely
        estimatesOutput <- apply(estimates[, -1], 2, round, decimal + 1)
        estimatesOutput <- data.frame(term = as.character(estimates$term), estimatesOutput, stringsAsFactors = FALSE)
        outputList$results2 <- estimatesOutput
    }
    
    if (showEffectSizesTable) {
        
        # get all other effect sizes
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes, stringsAsFactors = FALSE)
        
    }
    
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}

#### examples ####
# model1 <- lm(weight ~ Time, data = ChickWeight)
# model2 <- lm(-weight ~ Time, data = ChickWeight)
# model3 <- lm(mpg ~ drat, data = mtcars)
# # summary(model3)
# # reportLM(model1)
# reportLM(model1)
# # reportLM(model3, 2)
# # summary(model2)
# reportLM(model3)
# reportLM(model3, showTable = T, decimal = 3)
# reportLM(model3, showTable = T, showEffectSizesTable = T, decimal = 3)
# reportLM(model3, showTable = F, showEffectSizesTable = T, decimal = 2)
# summaryh(model3)