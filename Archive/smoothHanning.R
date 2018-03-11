smoothHanning <- function(x, windowLength = 11, window = 'hanning') {
    #prepare signal x by introducing reflected copies of signal with window size in both ends so that transient parts are minimized in the beginning and end part of the output signal
    signal <- c(rev(x[2:windowLength]), #
                x,
                rev(x[( length(x) - windowLength + 2):length(x)]))
    
    window <- signal::hanning(windowLength) #hanning window coefficients from signal package
    smoothedSignal <- stats::filter(signal, filt = window/sum(window))
    smoothedSignal <- smoothedSignal[-c(1:10, (length(smoothedSignal)-9):length(smoothedSignal) )] #remove begin and end padding
    if (length(x) == length(smoothedSignal)) {
        message('Original signal and smoothed signal equal length.')
    } else {
        warning("Smoothed signal is not of same length!")
    }
    return(smoothedSignal)
}
