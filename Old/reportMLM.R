reportMLM <- function(model, decimal = 2, showTable = FALSE, showEffectSizesTable = FALSE) {
    
    # ensure significant digits with sprintf
    digits <- paste0("%.", decimal, "f") # e.g, 0.10 not 0.1, 0.009, not 0.01
    pdigits <- paste0("%.", decimal + 1, "f") # p values always 1 more digit
    
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
        
        # effect size semi partial R (Edwards et al., 2008)
        anovaModel <- data.frame(anova(model))
        colnames(anovaModel) <- tolower(colnames(anovaModel))
        Fs <- anovaModel$f # F-values for each effect (marginal = type 3 SS with Satterthwaite (requires lmerTest package))
        numDF <- anovaModel$numdf #numerator DFs
        denDF <- anovaModel$dendf #denominator DFs
        semiPartialREffect <- (numDF / denDF * Fs) / (1 + (numDF / denDF * Fs)) # effect sizes
        
        if (length(semiPartialREffect) == nrow(estimates)) {
            estimates$es.partR2 <- semiPartialREffect
        } else if ( (nrow(estimates) - length(semiPartialREffect)) == 1 ) {
            estimates$es.partR2 <- c(NA, semiPartialREffect)
        }
        
        # piecewiseSEM (Nakagawa & Schielzeth, 2013)
        rsquareds <- sem.model.fits(model)
        estimates$es.R2marginal <- c(NA, rsquareds$Marginal)
        estimates$es.R2conditional <- c(NA, rsquareds$Conditional)
        
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
# library(lme4); library(lmerTest)
# model1 <- lmer(weight ~ Time  + (1 + Time | Chick), data = ChickWeight)
# model <- model1
# # model2 <- lmer(-weight ~ Time  + (1 + Time | Chick), data = ChickWeight)
# summary(model1)
# reportMLM(model1)
# reportMLM(model1, 2, T, T)
# reportMLM(model1, 2, F, T)
# reportMLM(model1, 2, T, F)
# 
# summary(model2)
# reportMLM(model2)

# library(nlme)
# model3 <- lme(weight ~ Time, random = ~ 1 + Time | Chick, data = ChickWeight)
# reportMLM(model3, 3, T, T)
# reportMLM(model3)
# reportMLM(model1)
# class(model3)[1]

# summaryh(model1)
# summaryh(model2)
# summaryh(model3)