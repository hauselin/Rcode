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
rcontrast <- function(t, df) {#
    return( sqrt(t^2 / (t^2 + df)) )
}

partialEta <- function(SSE, SSR) {
    return( SSE / (SSE + SSR) )
}



es <- function(d = NULL, r = NULL, R2 = NULL, f = NULL, oddsratio = NULL, logoddsratio = NULL, auc = NULL) {
    
    esList <- vector("list", 7)
    names(esList) <- c("d", "r", "R2", "f", "oddsratio", "logoddsratio", "auc")
    # auc calculations might be off...
    
    if (is.numeric(d)) {
        esList$d <- d
        esList$f <- esList$d / 2
        esList$r <- esList$d / sqrt(esList$d^2 + 4) # assumes equal sample size
        esList$R2 <- esList$r^ 2
        esList$oddsratio <- exp(esList$d / (sqrt(3) / pi))
        esList$logoddsratio <- esList$d / (sqrt(3) / pi)
        esList$auc <- pnorm(esList$d, 0, 1)
    } else if (is.numeric(r)) {
        esList$d <- (2 * r) / (sqrt(1 - r^2))
        esList$r <- r
        esList$f <- esList$d / 2
        esList$R2 <- esList$r^ 2
        esList$oddsratio <- exp(esList$d / (sqrt(3) / pi))
        esList$logoddsratio <- esList$d / (sqrt(3) / pi)
        esList$auc <- pnorm(esList$d, 0, 1)
    } else if (is.numeric(f)) {
        esList$d <- f * 2
        esList$r <- esList$d / sqrt(esList$d^2 + 4) # assumes equal sample size
        esList$f <- f
        esList$R2 <- esList$r^ 2
        esList$oddsratio <- exp(esList$d / (sqrt(3) / pi))
        esList$logoddsratio <- esList$d / (sqrt(3) / pi)
        esList$auc <- pnorm(esList$d, 0, 1)
    } else if (is.numeric(R2)) {
        esList$r <- sqrt(R2)
        esList$d <- (2 * esList$r) / (sqrt(1 - esList$r^2))
        esList$f <- esList$d / 2
        esList$R2 <- R2
        esList$oddsratio <- exp(esList$d / (sqrt(3) / pi))
        esList$logoddsratio <- esList$d / (sqrt(3) / pi)
        esList$auc <- pnorm(esList$d, 0, 1)
    }
    
    return(esList)
    
}

# es(d = 0.3)
# es(r = 0.3)
# es(f = 0.3)
# es(R2 = 0.6)
