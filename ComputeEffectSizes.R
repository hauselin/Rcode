omegaES <- function(SSM, SSR, dfM, dfR) {
        MSR <- SSR/dfR
        omegaSquared <- (SSM - dfM * MSR) / (SSM + SSR + MSR)
        #some people report omegaSquared, rather than omega
        #omega tends to be smaller than r because it adjusts for population estimation
        #.01 = small, .06 = medium, .14 = large (Kirk, 1996)
        omega <- sqrt(omegaSquared)
        #omega is simlar to r in terms of interpretation
        #.1 = small, .3 = medium, .5 = large
        return(omega)
}

t_r <- function(t, N, predictors) {
        return( sqrt(t^2 / (t^2 + N - predictors - 1)) )
}


# effect size for each contrast
rcontrast <- function(t, df) {
    return( sqrt(t^2 / (t^2 + df)) )
}

partialEta <- function(SSE, SSR) {
    return( SSE / (SSE + SSR) )
}