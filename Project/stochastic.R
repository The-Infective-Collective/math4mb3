library(Rcpp)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

nsim <- 100
R0vec <- seq(1, 20, by=0.2)
mvec <- c(0.001, 0.1)

reslist <- vector('list', length(mvec))

set.seed(101)
for (m in mvec) {
    
    M <- matrix(c(1-m, m, m, 1-m), 2, 2)
    
    subreslist <- vector('list', length(R0vec))
    
    for (R in R0vec) {
        print(paste(m,R, sep=","))
        
        pp <- base.params
        pp[["R0"]] <- R
        
        subsubreslist <- vector('list', nsim)
        
        for (i in 1:nsim) {
            
            init <- initfun(pp, 2, T)
            
            df <- SIRmodel_npatch_stochastic(pp, init, M, term_time)
            
            zero1 <- which(df$I[,1]==0)
            zero2 <- which(df$I[,2]==0)
            
            subsubreslist[[i]] <- data.frame(
                sim=i,
                R0=R,
                local=(any(df$I[,1]==0) || any(df$I[,2]==0)),
                global=(any(df$I[,1]==0 & df$I[,2]==0)),
                rescue1=length(which(diff(zero1) != 1)),
                rescue2=length(which(diff(zero2) != 1))
            )
            
        }
        
        subreslist[[which(R0vec==R)]] <- do.call("rbind", subsubreslist)
    }
    
    ss <- do.call("rbind", subreslist)
    ss$m <- m
    
    reslist[[which(mvec==m)]] <- ss
}

save("reslist", file="stochastic.rda")
