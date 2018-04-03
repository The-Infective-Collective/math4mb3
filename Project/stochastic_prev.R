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

