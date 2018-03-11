reportAOV <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {

    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    pdigits <- paste0("%.", decimal + 1, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01

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
    estimates$es.f <- cohens_f(model)
    estimates$es.r <- es(f = estimates$es.f, msg = F, decimal = decimal)$r

    # make a copy of estimates and convert to correct dp
    estimatesCopy <- estimates[, -1]
    estimatesRound <- estimatesCopy
    estimatesRound[estimatesCopy >= 0.01] <- round(estimatesRound[estimatesCopy >= 0.01], decimal)
    estimatesRound[estimatesCopy >= 0.01] <- sprintf(digits, estimatesCopy[estimatesCopy >= 0.01])
    estimatesRound[estimatesCopy < 0.01] <- sprintf(pdigits, estimatesCopy[estimatesCopy < 0.01])

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

    if (length(outputList) > 1) {
        return(outputList)
    } else {
        return(formattedOutputDf)
    }

}
