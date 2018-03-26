library(Rcpp)
sourceCpp("SIRmodel_onepatch.cpp")

N <- 1

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

pdf("bifurcation_R0.pdf", width=10, height=8)

plot(NA, 
     xlim=c(0, R0max), 
     ylim=c(1e-20, 1), 
     log="y",
     xlab="Reproductive number",
     ylab="Incidence")

for (R in R0) {
    params <- list(
        b0=17*365/5/N,
        b1=0.25,
        gamma=365/5,
        mu=0.02,
        a=2*pi,
        dt=1/3650,
        nsteps=3650*400
    )
    
    L <- lapply(initlist, function(x) {
        init <- x*N
        
        ## run the model once
        res <- SIRmodel(params, init)
        
        inc <- tail(res$I[res$time%%1==0], -300)/N
        
        points(rep(R, length(inc)), inc, pch=".")
    })
}

dev.off()

b1max <- 1
b1by <- 0.001

b1 <- seq(0, b1max, by=b1by)

pdf("bifurcation_amplitude.pdf", width=10, height=8)

plot(NA, 
     xlim=c(0, b1max), 
     ylim=c(1e-50, 1), 
     log="y",
     xlab="Seasonal amplitude",
     ylab="Incidence")

for (b in b1) {
    params <- list(
        b0=1250/N,
        b1=b,
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
        
        inc <- tail(res$incidence[res$time%%1==0], -80)/N
        
        points(rep(b, length(inc)), inc, pch=".")
    })
}

dev.off()
