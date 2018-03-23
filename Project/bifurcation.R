library(Rcpp)
sourceCpp("SIRmodel_onepatch.cpp")

N <- 10000

R0max <- 100
R0by <- 0.1

R0 <- seq(1, R0max, by=R0by)

initlist <- list(
    c(S=0.05, I=0.00001, R=1-0.05-0.00001),
    c(S=0.04, I=0.00001, R=1-0.04-0.00001),
    c(S=0.03, I=0.00001, R=1-0.03-0.00001),
    c(S=0.02, I=0.00001, R=1-0.02-0.00001),
    c(S=0.01, I=0.00001, R=1-0.01-0.00001)
)

pdf("bifurcation.pdf", width=8, height=6)

plot(NA, 
     xlim=c(0, R0max), 
     ylim=c(1e-20, 1), 
     log="y",
     xlab="Reproductive number",
     ylab="Incidence")

for (R in R0) {
    params <- list(
        b0=R*365/5/N,
        b1=0.25,
        gamma=365/5,
        mu=0.02,
        a=2*pi,
        dt=1/365,
        nsteps=365*100
    )
    
    L <- lapply(initlist, function(x) {
        init <- x*N
        
        ## run the model once
        res <- SIRmodel(params, init)
        
        inc <- tail(res$incidence[res$time%%1==0], -50)/N
        
        points(rep(R, length(inc)), inc, pch=".")
    })
}

dev.off()
