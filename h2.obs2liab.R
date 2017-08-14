## We need the inverse complementary error function to do the conversion
## Here's the function
erfcinv <- function(x) { qnorm(x/2, lower = FALSE)/sqrt(2) }

## Convert observed h2 to liability scale
## prevalence =  estimated trait prevalence
## cases = total number of cases, controls = total number of controls
## h2obs = h2 on the observed scale

h2liab <- function(prevalence, cases, controls, h2obs) {

       P <- cases/(cases + controls)
       Z <- sqrt(2)*erfcinv(2*prevalence)
       f <- 1/sqrt(2*pi) * exp(-Z^2/2)
       obs_to_liab <- (prevalence*(1-prevalence))^2/(P*(1-P))/f^2
       h2.liab <- h2obs * obs_to_liab

       return(h2.liab)

}
