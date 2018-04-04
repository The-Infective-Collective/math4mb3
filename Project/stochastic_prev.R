library(Rcpp)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

pp <- base.params
pp[["R0"]] <- 8.4

init <- initfun(pp, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df2 <- SIRmodel_npatch_stochastic(pp, init, M, term_time)

pdf("localext.pdf", width=8, height=6)

plot(df2$time, 
     df2$I[,2], 
     col=1, 
     type="l", 
     ylim=c(0, 8000), xlim=c(20, 40),
     xlab="Time (years)",
     ylab="Prevalence",
     main="Local Extinction Example")
lines(df2$time, df2$I[,1], col="blue")
points(df2$time[df2$I[,2]==0], rep(0, sum(df2$I[,2]==0)), col=2)
points(df2$time[df2$I[,1]==0], rep(0, sum(df2$I[,1]==0)), col=3)
legend(
    "topright",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

dev.off()

pp2 <- base.params
pp2[["R0"]] <- 2.0

init <- initfun(pp2, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df3 <- SIRmodel_npatch_stochastic(pp2, init, M, term_time)

pdf("globalext.pdf", width=8, height=6)

plot(df3$time, 
     df3$I[,2], 
     col=1, 
     type="l", 
     ylim=c(0, 6000), xlim=c(20, 50),
     xlab="Time (years)",
     ylab="Prevalence",
     main="Global Extinction Example")
lines(df3$time, df3$I[,1], col="blue")
points(df3$time[df3$I[,1]==0], rep(0, sum(df3$I[,1]==0)), col=3)
points(df3$time[df3$I[,2]==0], rep(0, sum(df3$I[,2]==0)), col=2)
legend(
    "topright",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

dev.off()

pp3 <- base.params
pp3[["R0"]] <- 17.0

init <- initfun(pp3, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df4 <- SIRmodel_npatch_stochastic(pp3, init, M, term_time)

pdf("asynchrony_lowm.pdf", width=8, height=6)

plot(df4$time, 
     df4$I[,2], 
     col=1, 
     type="l", 
     ylim=c(0, 6000), xlim=c(0, 50),
     xlab="Time (years)",
     ylab="Prevalence",
     main="Observed asynchrony with low m")
lines(df4$time, df4$I[,1], col="blue")
points(df4$time[df4$I[,2]==0], rep(0, sum(df4$I[,2]==0)), col=2)
points(df4$time[df4$I[,1]==0], rep(0, sum(df4$I[,1]==0)), col=3)
legend(
    "topright",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

dev.off()

pp4 <- base.params
pp4[["R0"]] <- 9.0

init <- initfun(pp, 2, T)
m <- 0.5
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df5 <- SIRmodel_npatch_stochastic(pp4, init, M, term_time)

#pdf("asynchrony_lowm.pdf", width=8, height=6)

plot(df5$time, 
     df5$I[,2], 
     col=1, 
     type="l", 
     ylim=c(0, 6000), xlim=c(0, 50),
     xlab="Time (years)",
     ylab="Prevalence",
     main="Observed asynchrony with low m")
lines(df5$time, df5$I[,1], col="blue")
points(df5$time[df5$I[,2]==0], rep(0, sum(df5$I[,2]==0)), col=2)
points(df5$time[df5$I[,1]==0], rep(0, sum(df5$I[,1]==0)), col=3)
legend(
    "topright",
    legend=c("Patch 1", "Patch 2"),
    col=c("blue", 1),
    lty=1
)

dev.off()

