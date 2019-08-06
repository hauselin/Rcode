# written by Hause Lin for Johnny Dubois

warning("You now have a new function in your environnment called shorts_in_public(). Type shorts_in_public() to begin.")

shorts_in_public <- function() {
    
    question <- "Should I wear shorts in public? Designed and made for Johnny Dubois."
    message(question)
    cat("Respond yes (y) or no (n) to the following questions to figure out.")

    r1 <- readline("Are you under 12? ")
    if (r1 == "y") {
        message("Have fun and wear your shorts, little man.")
    } else {
        r2 <- readline("Are you at the beach? ")
        if (r2 == "y") {
            message("Sure, go ahead and wear your shorts and enjoy your swim...")
        } else {
            r3 <- readline("Are you playing sports? ")
            if (r3 == "y") {
                message("Yea, OK, since you insist. Put on your shorts...")
            } else {
                r4 <- readline("Are you in North Africa helping the British defeat Erwin Rommel? ")
                if (r4 == "y") {
                    cat("Fine, weirdo... Rock your shorts.")
                } else {
                    warning("PUT SOME DAMN PANTS ON, JOHNNY!")
                }
            }
        }
    }

}
