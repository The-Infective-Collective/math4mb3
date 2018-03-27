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

df <- SIRmodel_npatch(base.params, initial, base.M, seasonal_cosine)

plot(df$I[,1], type="l")
lines(df$I[,2], col=2)

## or different values for seasonal amplitude

pp <- base.params

pp[["b1"]] <- 0.3

df2 <- SIRmodel_npatch(pp, initial, base.M, seasonal_cosine)

plot(df2$I[,1], type="l", xlim=c(0, 30000))
lines(df2$I[,2], col=2)


