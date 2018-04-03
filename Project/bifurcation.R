library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

L <- load("stochastic.rda")

nsim <- 20
R0vec <- seq(1, 20, by=0.2)

set.seed(101)
init <- data.frame(
    S=seq(from=0, to=0.1, length.out=nsim) * base.params[["pop"]],
    I=seq(from=0, to=0.0001, length.out=nsim) * base.params[["pop"]]
)

init[] <- apply(init, 2, sample)
init$R <- base.params[["pop"]] - (init$S + init$I)

blist <- vector('list', length(R0vec))

for (R in R0vec) {
    
    pp <- base.params
    pp[["R0"]] <- R
    pp[["nsteps"]] <- 365 * 4000 
    
    blist[[which(R0vec==R)]] <- lapply(1:nsim, function(x){
        print(paste(R, x))
        
        ii <- init[x,]
        
        df <- SIRmodel_npatch(pp, ii, matrix(1), term_time)
        
        prev <- tail(df$I[df$time%%1==0], 100)/pp[["pop"]]
        
        data.frame(
            prevalence=prev,
            R0=R
        )
    })
    
    save("blist", file="bifurcation.rda")
}

save("blist", file="bifurcation.rda")
