## How to use a new function

library(Rcpp)
source("params.R")

sourceCpp("SIRmodel_npatch.cpp")

## initial value takes vector of lists (need to specify initial value for each patch)
initial <- list(
    S = rep(base.init$S, 2) + rnorm(2),
    I = rep(base.init$I, 2) + rnorm(2),
    R = rep(base.init$R, 2) + rnorm(2)
)

df <- SIRmodel_npatch(base.params, initial, base.M, term_time)

plot(df$I[,1], type="l")
lines(df$I[,2], col=2)

## If you want to use a different parameter

pp <- base.params

pp[["b1"]] <- 0.2
df2 <- SIRmodel_npatch(pp, initial, base.M, term_time)

plot(df2$I[,1], type="l")
lines(df2$I[,2], col=2)

## or a different dispersal matrix

new.M <- matrix(c(0.999, 0.001, 0.001, 0.999), 2, 2)
initial3 <- list(
    S=c(0.05, 0.07) * 1e6,
    I=c(0.0001, 0.0001) * 1e6,
    R=(1-0.0001-c(0.05, 0.07)) * 1e6
)

df3 <- SIRmodel_npatch(base.params, initial3, new.M, term_time)

plot(df3$I[,1], type="l", xlim=c(30000, 35000))
lines(df3$I[,2], col=2)

