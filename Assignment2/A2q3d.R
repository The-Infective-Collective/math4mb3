## Load Philadelphia P&I data

datafile <- "pim_us_phila_city_1918_dy.csv"
philadata <- read.csv(datafile)
philadata$date <- as.Date(philadata$date)


## Define function that simulates mortality from prevalence curve
## (Prevalence curve is just I(t) solution curve)

mortality_simulation <- function(S0, I0, beta, gamma,
                                 eta, tau, tmax=nrow(philadata)) {
  soln <- as.data.frame(ode(y=c(x=S0, y=I0),
                            times = 1:(tmax+tau),
                            func = SIR.vector.field,
                            parms = c(beta=beta, gamma=gamma)))
  
  data.frame(date=philadata$date,
             pim=tail(soln$y, -tau)/eta)
}

## Trial and error with values

S0 <- 1-1e-7
I0 <- 1e-7
beta <- 0.55
gamma <- 0.25
eta <- 0.00025
tau <- 14

## Run simulation, plot solution curve vs. given data points

mfit <- mortality_simulation(S0=S0, I0=I0,
                             beta=beta, gamma=gamma,
                             eta=eta, tau=tau)
plot(mfit, type="l", lwd=3, col="blue", ylim=c(0, 800),
     xlab="Date", ylab="P&I Deaths")
points(philadata)

## Plot legend

legend("topright", legend=c("data", "fit"), lty=c(NA, 1),
       pch=c(1, NA), col=c("black", "blue"), lwd=2)

## Calculate optimal values for R_0 and mean infectious period

R_0 <- beta/gamma
mean_infectious_period <- 1/gamma
print(c(R_0, mean_infectious_period))

##So optimal fit occurs when:
    ## I0 = 1e-7, S0 = 1 - 1e-7
    ## R_0 = 2.2, mean infectious period = 4
    ## eta = 0.00025, tau = 14