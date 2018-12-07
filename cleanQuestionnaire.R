# Last modified by Hause Lin 24-10-17 6:18 PM
# Function to reverse-code questionnaires (wide form data).
# Computes overall and subscale (if any) mean, and reliability. 
# Function will first ask you to verify if you've specified your reverse-coded items and subscales correctly. If correctly specified, press 1 to continue. 
# Function returns a list containing wide form data, long form data, and reliability results.
# Parameters
# data: data frame
# subjectCol: column number indicating subject id (default = 1)
# scaleName: give your scale a name in characters (for example 'bigfive')
# scaleMin: scale minimum value (default = 1)
# scaleMax: scale maximu value (default = 7)
# subscales: provide subscale items in a list (see example below) (default = empty list/no subscales)
# itemsToReverse: items to reverse-code (default = nothing) 
# checkReliability: check reliability alpha (default = TRUE)

# WARNING: This function is still work in progress. If you input data is in long form, it might not work as intended!

# Example (BIS/BAS scale with subscales bis, basDrive, basFun, and basReward)
# scales$bisbas <- cleanQuestionnaire(data = bisbas, scaleMin = 1, scaleMax = 4, form = 'wide', scaleName = 'bisbas', subscales = list(bis = c(1, 6, 10, 13, 15, 18, 20), basDrive = c(2, 7, 9, 17), basFun = c(4, 8, 12, 16), basReward = c(3, 5, 11, 14, 19)), itemsToReverse = c(2:17, 19:20), checkReliability = T)


# install packages if necessary
packages <- c("tidyverse", "data.table", "dtplyr")
toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(toInstall)) install.packages(toInstall)
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE) 
library(tidyverse); library(data.table); library(dtplyr)
rm(packages); rm(toInstall)

# cleanQuestionnaire Function
cleanQuestionnaire <- function(data, subjectCol = 1, scaleName, scaleMin = 1, scaleMax = 7, subscales = list(), itemsToReverse = c(), checkReliability = TRUE) {
    
    message("Expecting wide data with columns in this order:\n1. participant number\n2, 3, 4... questions")
    
    data <- tbl_dt(data) # data.table
    participantVariable <- colnames(data)[subjectCol] # name of participant variable
    setnames(data, participantVariable, 'pNo')
    dataWithoutPNo <- select(data, -subjectCol) # data without participant variable
    dataWithoutPNo[, names(dataWithoutPNo) := lapply(.SD, as.character)] # in case it's factor, convert to character
    dataWithoutPNo[, names(dataWithoutPNo) := lapply(.SD, as.numeric)] # convert to numeric
    
    #### process wideform data ####
    
    scaleItems <- colnames(dataWithoutPNo) # variable names (scale item names)
    
    itemNamesToReverse <- c()
    if (length(itemsToReverse) > 0) { # if there are items to reverse
        itemNamesToReverse <- colnames(select(dataWithoutPNo, itemsToReverse)) # variable (item) names that need to be reverse coded
    } 
    
    # dataframe that indicates which items to reverse and subscales (for verification before processing)
    checkScaleItems <- tbl_dt(data_frame(scale = scaleName, itemName = colnames(dataWithoutPNo), item = as.numeric(1:ncol(dataWithoutPNo)), toReverse = 'no'))
    checkScaleItems[item %in% itemsToReverse, toReverse := 'yes']
    
    # add subscale name variable
    checkScaleItems$subscale <- scaleName
    if (length(subscales) > 0) { # if there are subscales, add them to checkScaleItems so we can verify before processing
        for (subscaleIdx in 1:length(subscales)) {
            # for each row/item, assign column 5 (subscale) the corresponding subscale name
            set(checkScaleItems, i = as.integer(subscales[[subscaleIdx]]), j = 5L, value = names(subscales)[subscaleIdx])
        }
    }

    # convert from wide to long form
    dataLong <- arrange(tbl_dt((gather(data, item, score, -subjectCol))), pNo)
    
    #### ask for confirmation to proceed ####
    
    checkScaleItems %>% arrange(subscale, item) %>% print(n = Inf)
    confirm <- ''
    # ask for confirmation (enter 1 to proceed, 0 to stop)
    confirm <- readline("Press 1 to proceed or 0 to stop: ")
    if (confirm == '0') {
        stop('Please check the dataframe you have provided to the data parameter!')
    } else if (confirm == '1') {
        message('Computing scale mean...')
    } else {
        stop('Please press 0 or 1.')
    }
    
    #### after confirming ####

    #### reverse code items ####
    # ensure scores are numeric; create new scoreR variable to store reverse-coded items later on
    dataLong[, `:=` (score = as.numeric(score), scoreR = as.numeric(score))]
    dataLong[item %in% itemNamesToReverse, scoreR := (scaleMax + 1 - score)] # reverse code selected items
    dataLong[, scaleName := scaleName]
    dataLong$subscale <- ''
    
    
    #### process differently depending on whether there are subscales ####
    if (length(subscales) > 0) { # if there are subscales
        
        subscaleNames <- names(subscales)
        # assign subscale name to each item
        for (subscaleI in 1:length(subscales)) {
            subscaleItems <- scaleItems[subscales[[subscaleI]]] # subscale item names for this subscale
            dataLong[item %in% subscaleItems, subscale := subscaleNames[subscaleI]]
        }
        
        # overall scale mean
        overallMean <- dataLong[, .(m = mean(scoreR, na.rm = T), 
                                    stdev = sd(scoreR, na.rm = T), # if sd is 0, problematic
                                    rge = range(scoreR, na.rm = T)[2] - range(scoreR, na.rm = T)[1], # if range is 0, problematic
                                    items = .N), 
                                by = .(pNo)]
        overallMean$subscale <- 'overall'
        
        # subscale means
        subscaleMean <- dataLong[, .(m = mean(scoreR, na.rm = T), 
                                     stdev = sd(scoreR, na.rm = T), 
                                     rge = range(scoreR, na.rm = T)[2] - range(scoreR, na.rm = T)[1], 
                                     items = .N), 
                                 by = .(pNo, subscale)]
        
        # combine overall mean and subscale mean dataframes
        summaryLong <- arrange(bind_rows(overallMean, subscaleMean), pNo)
        
        # summaryLong$subscale <- paste(scaleName, summaryLong$subscale, sep = '_')
        colnames(summaryLong) <- c(participantVariable, paste(scaleName, colnames(summaryLong)[-subjectCol], sep = '_')) # add _ to column name
        
        # convert to wide form
        summaryWide <- reshape(select(summaryLong, 1:2, 6), timevar = colnames(summaryLong)[6], idvar = colnames(summaryLong)[1], direction = "wide", sep = '_')
        
        # store wide and long forms in list (this function returns this list as output)
        scaleM <- list(wide = summaryWide, long = summaryLong)
        
        
    } else if (length(subscales) == 0) { # if no subscale
        
        summaryLong <- dataLong[, .(m = mean(scoreR, na.rm = T), 
                                    stdev = sd(scoreR, na.rm = T), 
                                    rge = range(scoreR, na.rm = T)[2] - range(scoreR, na.rm = T)[1], 
                                    items = .N), 
                                by = .(pNo)]
        summaryLong$subscale <- 'overall'
        colnames(summaryLong) <- c(participantVariable, paste(scaleName, colnames(summaryLong)[-subjectCol], sep = '_'))
        summaryWide <- reshape(select(summaryLong, 1:2, 6), timevar = colnames(summaryLong)[6], idvar = colnames(summaryLong)[1], direction = "wide", sep = '_')
        scaleM <- list(wide = summaryWide, long = summaryLong)
    }
    
    
    
    #### reliability (alpha) ####
    
    # if checkReliability == TRUE, load psych package to compute alpha (reliability of each subscale)
    if (checkReliability) {
        
        packages <- c("psych")
        toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
        if (length(toInstall)) {
            install.packages(toInstall)
        } else {
            library(psych)
        }
        
        reliability <- list() #empty list to score reliabilty results
        if (length(subscales) > 0) { # if there are subscales
            
            for (subscaleIdx in 1:length(subscales)) {
                subscaleName <- names(subscales)[subscaleIdx]
                tempData <- dataLong[subscale == subscaleName, c("pNo", 'item', 'scoreR'), with = FALSE] %>% # subset subscale items, select pNo, item, and scoreR columns
                    spread(item, scoreR) %>% # spread data to wide form
                    select(-1) %>% # remove participant number column 
                    as.data.frame() # psych::alpha function takes only dataframe
                
                reliability[[subscaleName]] <- psych::alpha(tempData, check.keys = TRUE) # store reliability results in list
                
            }
            
        } else { # if no subscales
            tempData <- dataLong[, c("pNo", 'item', 'scoreR'), with = FALSE] %>% # select pNo, item, and scoreR columns
                spread(item, scoreR) %>% # spread data to wide form
                select(-1) %>% # remove participant number column 
                as.data.frame() # psych::alpha function takes only dataframe
            
            # uses scaleName variable as list item name
            reliability[[scaleName]] <- psych::alpha(tempData, check.keys = TRUE) # store reliability results in list
        }
        
    } 
    
    # add reliability results to list that wil be returned by function
    scaleM$reliability <- reliability
    
    return(scaleM)
}

#### test function ####













