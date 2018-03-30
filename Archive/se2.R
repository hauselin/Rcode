#' Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#'    data: a data frame.
#'    measurevar: the name of a column that contains the variable to be summariezed
#'    groupvars: a vector containing names of columns that contain grouping variables
#'    na.rm: a boolean that indicates whether to ignore NA's
#'    conf.interval: the percent range of the confidence interval (default is 95%)
#'    http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

#' Last modified by Hause Lin 10-03-18 21:25 hauselin@gmail.com

# install packages if necessary
packages <- c("dplyr", "data.table")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) {
    install.packages(toInstall)
} else {
    library(dplyr); library(data.table);
}
rm(packages); rm(toInstall)

se2 <- function (data = NULL, measurevar, groupvars = NULL, na.rm = TRUE, conf.interval = 0.95, toNumeric = TRUE, sampleSize = NULL) {
    # between-subjects SE
    # convert to datatable
    data <- data.table(data)
    
    # function to compute N without NAs
    length2 <- function(x, na.rm = FALSE) {
        if (na.rm) 
            sum(!is.na(x))
        else length(x)
    }
    
    resultsList <- list() # empty list to store results
    
    for (i in 1:length(measurevar)) {
        
        if (sum( data.frame(data)[, measurevar[i]] %in% c(Inf, -Inf)) > 0) { # if measurvar contains Inf or -Inf, stop the script
            stop(paste0("\nInf or -Inf is in ", measurevar[i], " variable"))
        }
        
        # compute mean by group
        datac <- data[, .(unlist(lapply(.SD, length2, na.rm = na.rm)), 
                          unlist(lapply(.SD, mean, na.rm = na.rm)),
                          unlist(lapply(.SD, sd, na.rm = na.rm))),
                      by = groupvars, .SDcols = measurevar[i]]
        
        setnames(datac, c(groupvars, "N", measurevar[i], "sd")) # rename column names
        datac <- data.table(datac)
        setkeyv(datac, groupvars) # sort table
        
        if (!is.null(sampleSize)) {
            if (is.vector(sampleSize)) {
                datac$N <- sampleSize
            } else {
                datac$N <- NULL
                datac <- left_join(datac, sampleSize)
                datac <- data.table(datac)
            }
        }
        
        datac[, se := sd / sqrt(N)] # compute standard error
        
        ciMult <- qt(conf.interval / 2 + 0.5, unlist(datac$N) - 1)
        datac[, ci := se * ciMult]
        
        if (toNumeric) {
            # convert columns to numeric class if possible, else, leave as character
            oldwarning <- getOption("warn")
            options(warn = -1)
            for (j in 1:(ncol(datac)-4)) { # exclude last few columns (outcome, sd, se, ci)
                if (sum(is.na(as.numeric(as.character(datac[[j]])))) == 0) {
                    datac[[j]] <- as.numeric(datac[[j]])
                } else {
                    datac[[j]] <- as.character(datac[[j]])
                }
            }
            options(warn = oldwarning)    
        }
        
        datac <- select(datac, groupvars, N, everything())
        resultsList[[measurevar[i]]] <- datac
           
    }
    
    if (length(measurevar) == 1) {
        return(resultsList[[measurevar[1]]])
    } else {
        return(resultsList)
    }
    
}

#### examples ####
# se2(data = mtcars, measurevar = "disp", groupvars = c("cyl"))
# se2(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "am"))
# se2(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "vs"))



## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's

# norm data (this function will only be used by seWithin, and won't have to be called directly)
normWithin2 <- function (data = NULL, idvar, measurevar, betweenvars = NULL, na.rm = TRUE) {
    
    data <- data.table(data)
    setkeyv(data, idvar) # sort by idvar
    
    data.subjMean <- data[, .(unlist(lapply(.SD, mean, na.rm = na.rm))), by = c(idvar, betweenvars), .SDcols = measurevar] # compute mean for each subject
    setnames(data.subjMean, c(idvar, betweenvars,'subjMean'))
    dataNew <- left_join(data, data.subjMean)
    dataNew <- data.table(dataNew)
    setkeyv(dataNew, c(idvar, betweenvars)) # sort
    
    measureNormedVar <- paste0(measurevar, "Normed")
    
    # dataNew <- data.frame(dataNew)
    # dataNew[, measureNormedVar] <- dataNew[, measurevar] - unlist(data[, "subjMean"]) + mean(data[, measurevar], na.rm = na.rm)
    
    dataNew[, (measureNormedVar) := get(measurevar) - subjMean + mean(get(measurevar), na.rm = T)]
    dataNew$subjMean <- NULL
    
    return(data.frame(dataNew))
}

#### examples ####
# normWithin2(data = sleep, idvar = "ID", measurevar = "extra", betweenvars = "group") %>% arrange(ID)










## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
#'    showNormed: whether to show the normed version of the outcome variable

# within-subjects CI (normed and un-normed versions)
seWithin2 <- function (data = NULL, measurevar, betweenvars = NULL, withinvars = NULL, idvar = NULL, na.rm = TRUE, conf.interval = 0.95, showNormed = FALSE) {
    
    data <- data.frame(data) # convert to data.frame

    # Check if betweenvars and withinvars are factors
    factorvars <- sapply(data[, c(betweenvars, withinvars), drop = FALSE], 
                         FUN = is.factor)
    # Ensure that the betweenvars and withinvars are factors
    if (!all(factorvars)) {
        nonfactorvars <- names(factorvars)[!factorvars]
        message("Automatically converting the following non-factors to factors: ", 
                paste(nonfactorvars, collapse = ", "))
        data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
    }
    
    resultsList <- list() # empty list to store results
    
    for (i in 1:length(measurevar)) {
        
        # if measurvar contains Inf or -Inf, stop the script
        if (sum( data[, measurevar[i]] %in% c(Inf, -Inf)) > 0) { 
            stop(paste0("\nInf or -Inf is in ", measurevar[i], " variable"))
        }
        
        # Get the grand means from the un-normed data
        datac <- se2(data, measurevar[i], groupvars = c(betweenvars, withinvars), na.rm = na.rm, conf.interval = conf.interval, toNumeric = FALSE)
        
        # Drop all the unused columns (these will be calculated with normed data)
        datac$N <- NULL
        datac$sd <- NULL
        datac$se <- NULL
        datac$ci <- NULL
        
        # Norm each subject's data
        ndata <- normWithin2(data, idvar, measurevar[i], betweenvars, na.rm)
        
        # This is the name of the new column
        measurevar_n <- paste(measurevar[i], "Normed", sep = "")
        
        if (is.null(betweenvars)) {
            sampleSize <- length(unique(data[, idvar]))
        } else if (!is.null(betweenvars)) {
            sampleSize <- data.table(data)[, .N, by = c(betweenvars, idvar)][, .N, by = betweenvars]
        }
        
        # Collapse the normed data - now we can treat between and within vars the same
        ndatac <- se2(ndata, measurevar_n, groupvars = c(betweenvars, withinvars), na.rm = na.rm, conf.interval = conf.interval, toNumeric = FALSE, sampleSize = sampleSize)
        
        # Apply correction from Morey (2008) to the standard error and confidence interval
        # Get the product of the number of conditions of within-S variables
        # nWithinGroups <- prod(vapply(ndatac[,withinvars, drop = FALSE], FUN = function(x) length(levels(x)), FUN.VALUE = numeric(1)))
        nWithinGroups <- nrow(distinct(select(data, withinvars)))
        correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
        ndatacTbl <- data.table(ndatac)
        # Apply the correction factor
        ndatacTbl$sd <- ndatacTbl$sd * correctionFactor
        ndatacTbl$se <- ndatacTbl$se * correctionFactor
        ndatacTbl$ci <- ndatacTbl$ci * correctionFactor

        # Combine the un-normed means with the normed results
        merged <- left_join(datac, ndatacTbl)
        merged <- select(merged, c(betweenvars, withinvars), N, everything())
        merged <- mutate_if(data.frame(merged), is.factor, as.character) #if factor, convert to character
        merged <- data.table(merged)
        message("Factors have been converted to characters.")
        
        # convert columns to numeric class if possible, else, leave as character
        oldwarning <- getOption("warn")
        options(warn = -1)
        for (j in 1:(ncol(merged)-4)) { # exclude last few columns (outcome, sd, se, ci)
            if (sum(is.na(as.numeric(as.character(merged[[j]])))) == 0) {
                merged[[j]] <- as.numeric(merged[[j]])
            } else {
                merged[[j]] <- as.character(merged[[j]])
            }
        }
        options(warn = oldwarning)
        
        # whether to show normed version
        if (showNormed == FALSE) {
            # print(measurevar_n)
            merged[, (measurevar_n) := NULL]
        } 
        
        merged <- arrange_(merged, c(betweenvars, withinvars)) # arrange by between then within variables
        resultsList[[measurevar[i]]] <- merged
        message(cat("Confidence intervals: ", conf.interval, sep = ""))
        
    }
    
    if (length(measurevar) == 1) {
        # print(resultsList[[measurevar[1]]])
        return(resultsList[[measurevar[1]]])
    } else {
        # print(resultsList)
        return(resultsList)
    }
    
}

#### examples ####
# seWithin2(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick", showNormed = T)
# seWithin2(data = ChickWeight, measurevar = "weight", withinvars = "Time", idvar = "Chick")
# seWithin2(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick", showNormed = TRUE)
# library(ggplot2)
# seWithin2(data = diamonds, measurevar = c("carat"), betweenvars = "cut", withinvars = "color", idvar = "clarity")
# a <- seWithin2(data = diamonds, measurevar = c("carat", "depth", "z"), betweenvars = "cut", withinvars = "color", idvar = "clarity")
# a$carat
# a$depth
# a$z