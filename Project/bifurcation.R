library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

L <- load("stochastic.rda")

rinit <- function(param) {
    with(as.list(param), {
        ll <- list(
            S=runif(1, 0, 0.1) * pop,
            I=runif(1, 0, 0.0001) * pop
        )
        
        ll$R <- pop - ll$S - ll$I
        
        ll
    })
}

nsim <- 40
R0vec <- seq(1, 20, by=0.2)

blist <- vector('list', length(R0vec))

set.seed(101)
for (R in R0vec) {
    print(R)
    
    pp <- base.params
    pp[["R0"]] <- R
    pp[["nsteps"]] <- 365 * 200 
    
    blist[[which(R0vec==R)]] <- replicate(nsim, {
        init <- rinit(pp)
        
        df <- SIRmodel_npatch(pp, init, matrix(1), term_time)
        
        prev <- tail(df$I[df$time%%1==0], 20)/pp[["pop"]]
        
        data.frame(
            prevalence=prev,
            R0=R
        )
    }, simplify=FALSE)
    
    save("blist", file="bifurcation.rda")
}

save("blist", file="bifurcation.rda")
