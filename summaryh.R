# Last modified by Hause Lin 12-03-18 22:57 hauselin@gmail.com

# install really useful packages
packages <- c("dplyr", "data.table", "piecewiseSEM", "compute.es", "sjstats")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) {
    install.packages(toInstall)
} else {
    library(dplyr); library(data.table); library(piecewiseSEM); library(compute.es); library(sjstats)
}
rm(packages); rm(toInstall)

cat("r: .10 (small), .30 (medium), .50 (large) (Cohen, 1992)\n")
cat("d: 0.20 (small), 0.50 (medium), .80 (large) (Cohen, 1992)\n")
cat("R2: .02 (small), .13 (medium), .26 (large) (Cohen, 1992)\n")

summaryh <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    options(scipen = 999) # disable scientific notation
    if (class(model)[1] == 'lm') { # Last modified by Hause Lin 09-03-18 09:10 hauselin@gmail.com
        reportLM(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable) #
    } else if (class(model)[1] %in% c("glm", "glmerMod")) { # Last modified by Hause Lin 09-03-18 09:10 hauselin@gmail.com
        reportGLM(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable) #
    }  else if (class(model)[1] %in% c("lmerMod")) { # Last modified by Hause Lin 09-03-18 09:10 hauselin@gmail.com
        message("Please install/load lmerTest package and then refit your model!")
    } else if (class(model)[1] == "merModLmerTest") { # Last modified by Hause Lin 09-03-18 09:10 hauselin@gmail.com
        reportMLM(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (class(model)[1] == "aov") {
        reportAOV(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (class(model)[1] == "anova") {
        reportAOV(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (class(model)[1] == "lme") {
        reportMLM(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (grepl("t-test", model$method, ignore.case = T)) {
        reportTtest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (grepl("Pearson's product-moment correlation", model$method, ignore.case = T)) {
        reportCortestPearson(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else if (grepl("Kendall's rank correlation tau", model$method, ignore.case = T)) {
        reportCortest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable) 
    } else if (grepl("Spearman's rank correlation rho", model$method, ignore.case = T)) {
        reportCortest(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable) 
    } else if (grepl("Pearson's Chi-squared test", model$method, ignore.case = T)) {
        reportCHISQ(model = model, decimal = decimal, showTable = showTable, showEffectSizesTable = showEffectSizesTable)
    } else {
        # modelClass <- class(model)[1]
        message(paste0("Class not supported yet! Contact Hause Lin: hause.lin@mail.utoronto.ca"))
    }
    # options(scipen = scientificNotationDefault) # set to default scientific notation settings
}







reportLM <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
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
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
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
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}














reportAOV <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: F(3, 10) = 39, p < .001, r = 0.32
    if (class(model)[1] == "anova") {
        estimates <- data.frame(model) # get estimates and put in dataframe
    } else if (class(model)[1] == "aov") {
        estimates <- data.frame(anova(model)) # get estimates and put in dataframe    
    }
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$dfResid <- estimates["Residuals", "df"]  # get model degrees of freedom
    colnames(estimates) <- c('df', 'sum.sq', 'mean.sq', 'f.value', 'p.value', 'df.resid') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    estimates <- estimates[estimates$term != "Residuals", ]
    
    # effect sizes
    esCohensf <- cohens_f(model) # calculate Cohen's f
    if (is.data.frame(esCohensf)) { # sometimes output is a dataframe; if so, extract cohens.f variable
        estimates$es.f <- esCohensf$cohens.f
    } else {
        estimates$es.f <- esCohensf
    }
    
    estimates$es.r <- es(f = estimates$es.f, msg = F, decimal = decimal)$r
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    estimatesRound$df.resid <- round(estimates$df.resid)
    
    formattedOutput <- paste0("F(", estimatesRound$df, ", ", estimatesRound$df.resid, ")",
                              " = ", estimatesRound$f.value, 
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
        effectSizes <- es(r = round(as.numeric(estimates$es.omega), decimal + 1), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes, stringsAsFactors = FALSE)
        
    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}





reportCHISQ <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: r(30) = 0.82, p < .001
    estimates <- data.frame(df = model$parameter, 
                            statistic = model$statistic,
                            p.value = model$p.value)
    rownames(estimates) <- NULL
    
    # effect sizes
    computeES <- chies(chi.sq = estimates$statistic, n = sum(model$observed), verbose = F)
    estimates$es.r <- computeES$r
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
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
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}






reportCortest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: r = 0.82, p < .001
    estimates <- data.frame(estimate = model$estimate, 
                            p.value = model$p.value,
                            statistic = model$statistic)
    
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.d <- es(r = estimates$estimate, msg = F, decimal = decimal)$d # d
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    formattedOutput <- paste0("r = ", estimatesRound$estimate, 
                              ", p ", estimatesRound$p.value)
    
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
        effectSizes <- es(r = abs(round(as.numeric(estimates$estimate), decimal + 1)), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(model$data.name), as.list(effectSizes), stringsAsFactors = FALSE)
        
    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}








reportCortestPearson <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: r(30) = 0.82, p < .001
    
    estimates <- data.frame(df = model$parameter, 
                            estimate = model$estimate, 
                            estimateLow = model$conf.int[1],
                            estimateUpper = model$conf.int[2],
                            p.value = model$p.value,
                            statistic = model$statistic)
    
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.d <- (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    
    formattedOutput <- paste0("r(", estimatesRound$df, ")",
                              " = ", estimatesRound$estimate, 
                              ", p ", estimatesRound$p.value)
    
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
        effectSizes <- es(r = abs(round(as.numeric(estimates$estimate), decimal + 1)), decimal = decimal, msg = F)
        outputList$effectSizes <- data.frame(term = as.character(model$data.name), as.list(effectSizes), stringsAsFactors = FALSE)
        
    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}




reportGLM <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: b = −2.88, SE = 0.32, z(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- df.residual(model) # get model degrees of freedom
    colnames(estimates) <- c('estimate', 'std.error', 'statistic', 'p.value', 'df') # rename columns
    estimates <- data.frame(term = effectNames, estimates)
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.oddsratio <- exp(estimates$estimate)
    estimates$es.r <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$r # r
    estimates$es.d <- es(oddsratio = estimates$es.oddsratio, msg = F, decimal = decimal)$d # d
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
    # fix p values
    estimatesRound$p.value <- round(estimates$p.value, decimal + 2)
    estimatesRound$p.value <- ifelse(estimatesRound$p.value < .001, "< .001", paste0("= ", sprintf(pdigits, estimatesRound$p.value)))
    estimatesRound$p.value <- gsub("= 0.", "= .", estimatesRound$p.value)
    
    # leave df as integers
    estimatesRound$df <- round(estimates$df)
    
    formattedOutput <- paste0("b = ", estimatesRound$estimate, 
                              ", SE = ", estimatesRound$std.error, 
                              ", z(", estimatesRound$df, ")", 
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
        effectSizes <- es(r = round(as.numeric(estimates$es.r), decimal + 1), decimal = decimal)
        outputList$effectSizes <- data.frame(term = as.character(estimates$term), effectSizes, stringsAsFactors = FALSE)
        
    }
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}



reportMLM <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: b = −2.88, SE = 0.32, t(30) = −8.92, p < .001, r = .85
    estimates <- data.frame(coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    colnames(estimates) <- c('estimate', 'std.error', 'df', 'statistic', 'p.value') # rename columns
    estimates <- data.frame(term = effectNames, estimates, stringsAsFactors = FALSE)
    rownames(estimates) <- NULL
    
    # effect size r (Kashdan & Steger, 2006)
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
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
        
        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f.value # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        semiPartialREffect <- data.frame(term = rownames(anovaModel), es.partR2 = semiPartialREffect, stringsAsFactors = F)
        for (i in 1:nrow(semiPartialREffect)) {
            estimates$es.partR2[grepl(semiPartialREffect$term[i], estimates$term)] <- semiPartialREffect[i, "es.partR2"]
        }
        
        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        rsquareds <- sem.model.fits(model)
        estimates$es.R2marginal <- rsquareds$Marginal
        estimates$es.R2conditional <- rsquareds$Conditional
        
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
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}






reportTtest <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    if (decimal <= 2) {
        pdigits <- paste0("%.", 3, "f")
    } else {
        pdigits <- paste0("%.", decimal, "f")
    }
    
    # example output: t(30) = 5.82, p < .001
    estimates <- data.frame(df = model$parameter, 
                            statistic = model$statistic,
                            p.value = model$p.value)
    rownames(estimates) <- NULL
    
    # effect sizes
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[abs(estimatesCopy) >= 0.01] <- round(estimatesRound[abs(estimatesCopy) >= 0.01], decimal)
    estimatesRound[abs(estimatesCopy) >= 0.01] <- sprintf(digits, estimatesCopy[abs(estimatesCopy) >= 0.01])
    estimatesRound[abs(estimatesCopy) < 0.01] <- signif(estimatesCopy[abs(estimatesCopy) < 0.01], digits = 1)
    # estimatesRound[abs(estimatesCopy) < 0.01] <- sprintf(pdigits, estimatesCopy[abs(estimatesCopy) < 0.01])
    
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
    options(scipen = 0) # enable scientific notation
    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }
    
}






es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL, decimal = 3, msg = TRUE) {
    
    # Last modified by Hause Lin 22-11-17 09:51 hauselin@gmail.com
    # effectsizes <- vector("list", 7) # list version
    effectsizes <- data.frame(matrix(NA, nrow = length(c(d, r, R2, f, oddsratio, logoddsratio, auc)), ncol = 7)) # dataframe version
    names(effectsizes) <- c("d", "r", "R2", "f", "oddsratio", "logoddsratio", "auc")
    # auc calculations might be off...
    
    if (length(c(d, r, R2, f, oddsratio, logoddsratio, auc)) < 1) { 
        stop("Please specify one effect size!")
    }
    
    if (is.numeric(d)) {
        if (msg) {message(paste0("d: ", d, " ")) }
        effectsizes$d <- d
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(r)) {
        if (msg) {message(paste0("r: ", r, " ")) }
        effectsizes$d <- (2 * r) / (sqrt(1 - r^2))
        effectsizes$r <- r
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(f)) {
        if (msg) {message(paste0("f: ", f, " ")) }
        effectsizes$d <- f * 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$f <- f
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(R2)) {
        if (msg) {message(paste0("R2: ", R2, " ")) }
        effectsizes$r <- sqrt(R2)
        effectsizes$d <- (2 * effectsizes$r) / (sqrt(1 - effectsizes$r^2))
        effectsizes$f <- effectsizes$d / 2
        effectsizes$R2 <- R2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(oddsratio)) {
        if (msg) {message(paste0("odds ratio: ", oddsratio, " "))}
        effectsizes$d <- log(oddsratio) * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- oddsratio
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(logoddsratio)) {
        if (msg) {message(paste0("log odds ratio: ", logoddsratio, " ")) }
        effectsizes$logoddsratio <- logoddsratio
        effectsizes$d <- effectsizes$logoddsratio * (sqrt(3) / pi)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$auc <- pnorm(effectsizes$d, 0, 1)
    } else if (is.numeric(auc)) {
        if (msg) {message(paste0("auc: ", auc, " ")) }
        effectsizes$auc <- auc
        effectsizes$d <- qnorm(auc, 0, 1)
        effectsizes$f <- effectsizes$d / 2
        effectsizes$r <- effectsizes$d / sqrt(effectsizes$d^2 + 4) # assumes equal sample size
        effectsizes$R2 <- effectsizes$r^ 2
        effectsizes$oddsratio <- exp(effectsizes$d / (sqrt(3) / pi))
        effectsizes$logoddsratio <- effectsizes$d / (sqrt(3) / pi)
    }
    
    return(round(effectsizes, decimal))
    
}


# es(d = 0.3)
# es(d = 0.3, r = 0.2)
# es(d = c(0.3, 0.4), r = 0.2, f = 0.5)
# es(d = c(0.2, 0.3, 0.4))
# es(d = c(0.2, 0.3))$r
# es(r = 0.5)
# es(f = 0.24)
# es(R2 = 0.6)
# es(oddsratio = 1.6)
# es(logoddsratio = 1.6)
# es(auc = .99)
# es()