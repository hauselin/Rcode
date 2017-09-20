reportLM <- function(model, decimal = 3, intercept = FALSE, format = 'apa', showTable = TRUE) {
    # APA format
    message("Summary statistics for linear models fitted with lm()")
    cat("r: .10 (small), .30 (medium), .50 (large) (Cohen, 1992)\nd: 0.20 (small), 0.50 (medium), .80 (large) (Cohen, 1992)\nR2: .02 (small), .13 (medium), .26 (large) (Cohen, 1992)")
    
    estimates <- data.frame(coef(summary(model))) # get estimates and put in dataframe
    effectNames <- rownames(estimates) # get names of effects
    colnames(estimates) <- tolower(colnames(estimates))
    estimates$df <- df.residual(model) # get model degrees of freedom
    estimates <- estimates[, c(1, 2, 5, 3, 4)] # sort columns
    colnames(estimates) <- c('estimate', 'se', 'df', 'statistic', 'p') # rename columns
    
    # effect sizes
    estimates$es.r <-  sqrt((estimates$statistic ^ 2 / (estimates$statistic ^ 2 + estimates$df))) # r
    estimates$es.d <-  (2 * estimates$statistic) / sqrt(estimates$df) # d
    
    # R2
    estimates$r.squared <- summary(model)$r.squared
    estimates$adj.r.squared <- summary(model)$adj.r.squared
    
    # print formatted results
    if (intercept) {
        startingRow <- 1 # if intercept == TRUE, then report intercept too
    } else {
        startingRow <- 2 # if intercept == FALSE, then don't report intercept
    }
    
    for (i in startingRow:nrow(estimates)) { # for each fixed effect, format output according to APA style
        
        # i <- i
        estimateV <- estimates[i, 'estimate']
        estimateV <- ifelse(abs(estimateV) < 0.01, round(estimateV, 3), round(estimateV, 2))
        estimateSign <- sign(estimateV) # sign (negative or positive)
        estimateV <- ifelse(abs(estimateV) < 0.01, 
                            sprintf('%.3f', abs(estimateV)), 
                            sprintf('%.2f', abs(estimateV))) # character (absolute value)
        if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative signs are dashes, not hyphens
            estimateV <- ifelse(estimateSign == -1, paste0("–", estimateV), estimateV) # if negative, add minus sign to character vector
        } else {
            estimateV <- ifelse(estimateSign == -1, paste0("-", estimateV), estimateV) # if negative, add minus sign to character vector
        }
        
        seV <- estimates[i, 'se']
        seV <- ifelse(abs(seV) < 0.01, round(seV, 3), round(seV, 2))
        seV <- ifelse(abs(seV) < 0.01, 
                      sprintf('%.3f', abs(seV)), 
                      sprintf('%.2f', abs(seV))) # character
        dfV <- round(estimates[i, 'df'])
        statisticV <- estimates[i, 'statistic']
        statisticSign <- sign(statisticV)
        statisticV <- ifelse(abs(statisticV) < 0.01, 
                             sprintf('%.3f', abs(statisticV)), 
                             sprintf('%.2f', abs(statisticV))) # character
        if (.Platform$OS.type == 'unix') { # if linux/mac, ensure negative signs are dashes, not hyphens
            statisticV <- ifelse(statisticSign == -1, paste0("–", statisticV), statisticV) # character
        } else {
            statisticV <- ifelse(statisticSign == -1, paste0("-", statisticV), statisticV) # character
        }
        
        pV <- estimates[i, 'p']
        rV <- estimates[i, 'es.r']
        
        if (estimates[i, 'p'] >= 0.001) { # if p value is >= .001, always print exact p value
            
            if (format == 'nn') { # if nature neuro format
                
                message(sprintf("%.0f: b = %s, SE = %s, t%.0f = %s, P = %.3f, r = %.2f",
                                i, estimateV, seV, dfV, statisticV, pV, rV))  
                
            } else { # else always return APA format
                
                message(sprintf("%.0f: b = %s, SE = %s, t(%.0f) = %s, p = %s, r = %s",
                                i, estimateV, seV, dfV, statisticV, substring(sprintf('%.3f', pV), 2), substring(sprintf('%.2f', rV), 2)))
                
            }
            
        } else if (estimates[i, 'p'] < 0.001) { # if p value < .001, then always print p < .001
            
            pV <- 0.001
            
            if (format == 'nn') {
                
                message(sprintf("%.0f: b = %s, SE = %s, t%.0f = %s, P < %.3f, r = %.2f",
                                i, estimateV, seV, dfV, statisticV, pV, rV))  
                
            } else {
                
                message(sprintf("%.0f: b = %s, SE = %s, t(%.0f) = %s, p < %s, r = %s",
                                i, estimateV, seV, dfV, statisticV, substring(sprintf('%.3f', pV), 2), substring(sprintf('%.2f', rV), 2)))
                
            }
            
        }
    }

    rownames(estimates) <- paste(1:nrow(estimates), rownames(estimates)) # add row number to summary table
    
    if (showTable) {
        return(round(estimates, decimal))    
    }
    
}

#### test function #####
# model1 <- lm(weight ~ Time, data = ChickWeight)
# model2 <- lm(-weight ~ Time, data = ChickWeight)
# summary(model1)
# reportLM(model1)

# summary(model2)
# reportLM(model2)