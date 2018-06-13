# cleanQuestionnairePsychopy Function
# Last modified by Hause Lin 24-10-17 6:19 PM
cleanQuestionnairePsychopy <- function(data, subjectCol = 1, scaleName, scaleMin = 1, scaleMax = 7, subscales = list(), itemsToReverse = c(), checkReliability = TRUE) {
    
    message("Expecting long/tidy data with four columns in this order:\n1. participant number\n2. item number in scale\n3. subscale name\n4. rating\n")
    
    # install packages if necessary
    packages <- c("dplyr", "data.table", "dtplyr", "tidyr")
    toInstall <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(toInstall)) install.packages(toInstall)
    library(tidyr); library(dplyr); library(data.table); library(dtplyr)
    rm(packages); rm(toInstall)
    
    data <- tbl_dt(data) # data.table
    participantVariable <- colnames(data)[subjectCol] # name of participant variable
    setnames(data, c('pNo', 'item', 'subscale', 'score'))
    data[, pNo := as.numeric(as.character(pNo))] # ensure participant number column is numeric
    data[, score := as.numeric(as.character(score))] # ensure rating column is numeric
    
    # determine how many questions/items in scale
    u <- data[, .N, by = pNo][, N] # number of items per participant 
    scaleN <- u[which.max(tabulate(match(u, unique(u))))] # modal number of items (in case participants didn't complete scales)
    # if (length(scaleN) != 1) {
    #     print(data[, .N, by = pNo], n = Inf)
    #     stop("Participants do not have the same number of items!")
    # }
    scaleItems <- data[1:scaleN, as.character(item)]
    subscalesNames <- data[1:scaleN, subscale]
    
    # determine names of items to reverse
    itemNamesToReverse <- c()
    if (length(itemsToReverse) > 0) {
        itemNamesToReverse <- scaleItems[itemsToReverse]
    } 
    
    checkScaleItems <- tbl_dt(data_frame(
        scale = scaleName,
        itemName = scaleItems, 
        item = as.numeric(1:length(scaleItems)),
        toReverse = 'no',
        subscale = subscalesNames))
    checkScaleItems[item %in% itemsToReverse, toReverse := 'yes']
    
    
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
    data[, `:=` (score = as.numeric(score), scoreR = as.numeric(score))]
    data[item %in% itemNamesToReverse, scoreR := (scaleMax + 1 - score)] # reverse code selected items
    data[, scaleName := scaleName]
    data$subscale <- ''
    
    #### process differently depending on whether there are subscales ####
    if (length(subscales) > 0) { # if there are subscales
        
        subscaleNames <- names(subscales)
        # assign subscale name to each item
        for (subscaleI in 1:length(subscales)) {
            subscaleItems <- scaleItems[subscales[[subscaleI]]] # subscale item names for this subscale
            data[item %in% subscaleItems, subscale := subscaleNames[subscaleI]]
        }
        
        # overall scale mean
        overallMean <- data[, .(m = mean(scoreR, na.rm = T), 
                                    stdev = sd(scoreR, na.rm = T), # if sd is 0, problematic
                                    rge = range(scoreR, na.rm = T)[2] - range(scoreR, na.rm = T)[1], # if range is 0, problematic
                                    items = .N), 
                                by = .(pNo)]
        overallMean$subscale <- 'overall'
        
        # subscale means
        subscaleMean <- data[, .(m = mean(scoreR, na.rm = T), 
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
        
        summaryLong <- data[, .(m = mean(scoreR, na.rm = T), 
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
                tempData <- data[subscale == subscaleName, c(participantVariable, 'item', 'scoreR'), with = FALSE] %>% # subset subscale items, select pNo, item, and scoreR columns
                    spread(item, scoreR) %>% # spread data to wide form
                    select(-1) %>% # remove participant number column 
                    as.data.frame() # psych::alpha function takes only dataframe
                
                reliability[[subscaleName]] <- psych::alpha(tempData, check.keys = TRUE) # store reliability results in list
                
            }
            
        } else { # if no subscales
            tempData <- data[, c(participantVariable, 'item', 'scoreR'), with = FALSE] %>% 
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
# cleanQuestionnairePsychopy(data, subjectCol = 1, scaleName, scaleMin = 1, scaleMax = 7, subscales = list(), itemsToReverse = c(), checkReliability = TRUE)










