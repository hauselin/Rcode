## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
se <- function (data = NULL, measurevar, groupvars = NULL, na.rm = TRUE, conf.interval = 0.95) {
    
    # install packages if necessary
    packages <- c("dplyr", "data.table", "dtplyr")
    toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(toInstall)) {
        install.packages(toInstall)
    } else {
        library(dplyr); library(data.table); library(dtplyr)
    }
    
    # convert to datatable and tibble
    data <- tbl_dt(data)
    
    # function to compute N without NAs
    length2 <- function(x, na.rm = FALSE) {
        if (na.rm) 
            sum(!is.na(x))
        else length(x)
    }
    
    resultsList <- list() # empty list to store results
    
    for (i in 1:length(measurevar)) {
        
        # compute mean by group
        datac <- data[, .(unlist(lapply(.SD, length2, na.rm = na.rm)), 
                          unlist(lapply(.SD, mean, na.rm = na.rm)),
                          unlist(lapply(.SD, sd, na.rm = na.rm))),
                      by = groupvars, .SDcols = measurevar[i]]
        
        setnames(datac, c(groupvars, "N", measurevar[i], "sd")) # rename column names
        
        setkeyv(datac, groupvars) # sort table
        
        datac[, se := sd / sqrt(N)] # compute standard error
        
        ciMult <- qt(conf.interval / 2 + 0.5, unlist(datac$N) - 1)
        datac[, ci := se * ciMult]
        
        resultsList[[measurevar[i]]] <- tbl_df(datac)
        
    }
    
    if (length(measurevar) == 1) {
        return(resultsList[[measurevar[1]]])
    } else {
        return(resultsList)
    }
    
}

#### test function ####
# se(data = mtcars, measurevar = "disp", groupvars = c("cyl"))
# se(data = mtcars, measurevar = c("mpg", "disp"), groupvars = c("cyl", "vs"))







## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's

# norm data (this function will only be used by seWithin, and won't have to be called directly)
normWithin <- function (data = NULL, idvar, measurevar, betweenvars = NULL, na.rm = TRUE) {
    
    # install packages if necessary
    packages <- c("dplyr", "data.table", "dtplyr")
    toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(toInstall)) {
        install.packages(toInstall)
    } else {
        library(dplyr); library(data.table); library(dtplyr)
    }
    
    data <- tbl_dt(data)
    setkeyv(data, idvar) # sort by idvar
    
    data.subjMean <- data[, .(unlist(lapply(.SD, mean, na.rm = na.rm))), by = c(idvar, betweenvars), .SDcols = measurevar] # compute mean for each subject
    setnames(data.subjMean, c(idvar, betweenvars,'subjMean'))
    dataNew <- left_join(data, data.subjMean)
    setkeyv(dataNew, c(idvar, betweenvars)) # sort
    
    measureNormedVar <- paste0(measurevar, "Normed")
    
    # dataNew <- data.frame(dataNew)
    # dataNew[, measureNormedVar] <- dataNew[, measurevar] - unlist(data[, "subjMean"]) + mean(data[, measurevar], na.rm = na.rm)
    
    dataNew[, (measureNormedVar) := get(measurevar) - subjMean + mean(get(measurevar), na.rm = T)]
    
    dataNew$subjMean <- NULL
    
    return(data.frame(dataNew))
}

#### test function ####
# data <- ChickWeight
# dvar <- "Chick"
# measurevar <- "weight"
# betweenvars <- "Diet"
# subjMean <- "subjMean"
# na.rm = T
# tbl_dt(normWithin(data = ChickWeight, idvar = "Chick", measurevar = "weight", betweenvars = "Diet"))
# normWithin(data = sleep, idvar = "ID", measurevar = "extra", betweenvars = "group") %>% arrange(ID)










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

# within-subjects CI (normed and un-normed versions)
seWithin <- function (data = NULL, measurevar, betweenvars = NULL, withinvars = NULL, 
                              idvar = NULL, na.rm = TRUE, conf.interval = 0.95) {
    
    # install packages if necessary
    packages <- c("dplyr", "data.table", "dtplyr")
    toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(toInstall)) {
        install.packages(toInstall)
    } else {
        library(dplyr); library(data.table); library(dtplyr)
    }
    
    data <- data.frame(data) # convert to data.frame
    
    # Check if betwenvars and withinvars are factors
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
        
        # Get the means from the un-normed data
        datac <- summarySE2(data, measurevar[i], groupvars = c(betweenvars, withinvars),
                            na.rm = na.rm, conf.interval = conf.interval)
        
        # Drop all the unused columns (these will be calculated with normed data)
        datac$sd <- NULL
        datac$se <- NULL
        datac$ci <- NULL
        
        # Norm each subject's data
        ndata <- normWithin(data, idvar, measurevar[i], betweenvars, na.rm)
        
        # This is the name of the new column
        measurevar_n <- paste(measurevar[i], "Normed", sep = "")
        
        # Collapse the normed data - now we can treat between and within vars the same
        ndatac <- se(ndata, measurevar_n, groupvars = c(betweenvars, withinvars),
                     na.rm = na.rm, conf.interval = conf.interval)
        
        # Apply correction from Morey (2008) to the standard error and confidence interval
        # Get the product of the number of conditions of within-S variables
        nWithinGroups <- prod(vapply(ndatac[,withinvars, drop = FALSE], FUN = function(x) length(levels(x)), FUN.VALUE = numeric(1)))
        correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
        
        ndatacTbl <- tbl_dt(ndatac)
        
        # Apply the correction factor
        ndatacTbl[, `:=` (sd = sd * correctionFactor, se = se * correctionFactor, ci = ci * correctionFactor)]
        
        # Combine the un-normed means with the normed results
        merged <- left_join(datac, ndatacTbl)
        merged <- mutate_if(merged, is.factor, as.character) #if factor, convert to character
        merged[order( unlist((merged[, 1])), decreasing =  F), ] #arrange by first column
        message("Factors have been converted to characters.")
        
        resultsList[[measurevar[i]]] <- merged
        
    }
    
    
    
    if (length(measurevar) == 1) {
        return(resultsList[[measurevar[1]]])
    } else {
        return(resultsList)
    }
    
}

#### test function ####
# seWithin(data = ChickWeight, measurevar = "weight", betweenvars = "Diet", withinvars = "Time", idvar = "Chick")
# library(ggplot2)
# seWithin(data = diamonds, measurevar = c("carat"), betweenvars = "cut", withinvars = "color", idvar = "clarity")
# a <- seWithin(data = diamonds, measurevar = c("carat", "depth", "z"), betweenvars = "cut", withinvars = "color", idvar = "clarity")
# a$carat
# a$depth
# a$z