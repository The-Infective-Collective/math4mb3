library(Rcpp)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

pp <- base.params
pp[["R0"]] <- 8.4

init <- initfun(pp, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(2)
df2 <- SIRmodel_npatch_stochastic(base.params, init, M, term_time)

#pdf("stochastic1.pdf", width=8, height=6)

plot(df2$time, 
     df2$I[,2], 
     col=1, 
     type="l", 
     ylim=c(0, 30000), xlim=c(0, 1),
     xlab="Time (years)",
     ylab="Prevalence")
lines(df2$time, df2$I[,1], col="blue")
points(df2$time[df2$I[,2]==0], rep(0, sum(df2$I[,2]==0)), col=2)
points(df2$time[df2$I[,1]==0], rep(0, sum(df2$I[,1]==0)), col=3)
legend(
    "topright",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

## STOCHASTIC DRIFT

pp2 <- base.params
pop <- pp2[["pop"]] <- 7e5

init <- list(
    S=c(0.05, 0.07) * pop,
    I=c(0.0001, 0.0001) *  pop,
    R=(1-0.0001-c(0.05, 0.07)) * pop
)

M <- matrix(c(0.999, 0.001, 0.001, 0.999), 2, 2)

set.seed(4)
df2 <- SIRmodel_npatch_stochastic(base.params, init, M, seasonal_cosine)

pdf("stochastic2.pdf", width=8, height=6)

plot(df2$time, 
     df2$I[,2], 
     col=1, 
     type="l", 
     xlab="Time (years)",
     ylab="Prevalence",
     xlim=c(50, 100))
lines(df2$time, df2$I[,1], col="blue")
points(df2$time[df2$I[,2]==0], rep(0, sum(df2$I[,2]==0)), col=2)
legend(
    "topleft",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

dev.off()

