# Last modified by Hause Lin 20-09-17 17:21

#### z-score outlier detection ####
outliersZ <- function(x, zCutOff = 1.96, replaceOutliersWith = NA, outlierIndices = FALSE, showZValues = FALSE, digits = 2) {
    # clean and identify outliers using Z-score cut-off method
    # zCutOff: value to use (1.96 is standard)
    # replaceOutliersWith: if value is an outlier, what to replace it with? NA by default
    # outlierIndices: return index of outlier
    # showZValues: if TRUE, will show z score of each value
    # digits: rounding digits
    
    # compute standard deviation (sample version n = n [not n-1])
    stdev <- sqrt(sum((x - mean(x, na.rm = T))^2, na.rm = T) / sum(!is.na(x)))
    # compute absolute Z values for each value
    absZ <- abs(x - mean(x, na.rm = T)) / stdev
    # subset data that has absZ greater than the zCutOff and replace them with replace
    # can also replace with other values (such as max/mean of data)
    x[absZ > zCutOff] <- replaceOutliersWith
    outliers <- length(x[absZ > zCutOff])
    
    if (showZValues) {
        message("Showing absolute z-scores for each value.")
        message(paste0(outliers, " outliers detected."))
        return(round(absZ, digits)) # if values == TRUE, return z score for each value
    } else if (outlierIndices) {
        message("Showing indices of outliers.")
        return(which(is.na(x)))
    } else {
        message(paste0(outliers, " outliers detected."))
        message(paste0("Outliers replaced with ", replaceOutliersWith))
        return(round(x, digits)) # otherwise, return values with outliers replaced
    }
}

#### median absolute deviation outlier detection ####
outliersMAD <- function(x, MADCutOff = 2.5, replaceOutliersWith = NA, showMADValues = FALSE, outlierIndices = FALSE, bConstant = 1.4826, digits = 2) {
    # clean and identify outliers using MAD cut-off method (see Leys et al., 2013)
    # x: vector
    # MADCutOFF: value to use (2.5 is recommended)
    # replaceOutliersWith: if value is an outlier, what to replace it with? NA by default
    # showMADValues:if TRUE, will show MAD for each value instead
    # outlierIndices: return index of outlier
    # bConstant: usually, b = 1.4826, a constant linked to the assumption of normality of the data, disregarding the abnormality induced by out- liers (Rousseeuw & Croux, 1993).
    # digits: rounding digits
    
    # compute number of absolute MADs away for each value: formula: abs( ( x - median(x) ) )/ mad(x)
    absMADAway <- abs((x - median(x, na.rm = T))/mad(x, constant = bConstant, na.rm = T))
    # subset data that has absMADAway greater than the MADCutOff and replace them with replace
    x[absMADAway > MADCutOff] <- replaceOutliersWith
    outliers <- length(x[absMADAway > MADCutOff])
    if (showMADValues) { # if values == TRUE, return number of mads for each value
        message("Showing absolute MAD from median for each value.")
        message(paste0(outliers, " outliers detected."))
        return(round(absMADAway, digits)) 
    } else if (outlierIndices) {
        message("Showing indices of outliers.")
        return(which(is.na(x)))
    } else {
        message(paste0(outliers, " outliers detected."))
        message(paste0("Outliers replaced with ", replaceOutliersWith))
        return(round(x, digits)) # otherwise, return original with outliers replaced
    }
}


#### test functions above ####
# x <- c(1, 3, 3, 6, 8, 10, 10, 1000)
# outliersZ(x, showZValues = F, zCutOff = 1, replaceOutliersWith = NA, outlierIndices = T)
# outliersMAD(x, showMADValues = F, replaceOutliersWith = NA, outlierIndices = T)
