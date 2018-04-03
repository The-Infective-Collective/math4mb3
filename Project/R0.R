library(Rcpp)
library(tidyverse)
source("params.R")
source("useful-fun.R")
sourceCpp("SIRmodel_npatch.cpp")


nsim <- 100
R0vec <- seq(1, 20, by=0.2)
mvec <- c(0.001, 0.01, 0.1, 0.5)

reslist <- vector('list', length(mvec))

set.seed(101)

for (m in mvec) {
  
  M <- matrix(c(1-m, m, m, 1-m), 2, 2)
  
  subreslist <- vector('list', length(R0vec))
  
  for (R in R0vec) {
    print(paste(m,R, sep=","))
    
    pp <- base.params
    pp[["R0"]] <- R
    
    subsubreslist <- vector("list", nsim)
    
 
    for (i in 1:nsim) {
      
      init <- r.init.science(2)
      
      df <- run_SIR(pp, init, M, term_time)

      incoherence2 <- apply(df[(nrow(df) - 2*365):nrow(df), c("I1", "I2")], 1, coherence_calc) %>% mean()
      incoherence5 <- apply(df[(nrow(df) - 5*365):nrow(df), c("I1", "I2")], 1, coherence_calc) %>% mean()
      
      subsubreslist[[i]] <- data_frame(
          S01 = init$S[1],
          S02 = init$S[2],
          I01 = init$I[1],
          I02 = init$I[2],
          R0=R,
          incoherence2 = incoherence2,
          incoherence5 = incorherence5
      )
    }
    
    subreslist[[which(R0vec ==R)]] <-do.call("rbind", subsubreslist)
  }
  
  ss <- do.call("rbind", subreslist)
  ss$m <- m
  
  reslist[[which(mvec==m)]] <- ss
}

save("reslist", file="stochastic.rda")
