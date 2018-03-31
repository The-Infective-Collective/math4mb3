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
            S=runif(1, 0, 1) * pop,
            I=runif(1, 0, 0.0001) * pop
        )
        
        ll$R <- pop - ll$S - ll$I
        
        ll
    })
}

base.params <- list(
    R0=17,
    pop=1e6,
    b1=0.25,
    gamma=365/13,
    mu=0.02,
    dt=1/365,
    nsteps=365*100
)

nsim <- 1
R0vec <- seq(1, 40, by=0.2)

blist <- vector('list', length(R0vec))

for (R in R0vec) {
    print(R)
    
    pp <- base.params
    pp[["R0"]] <- R
    pp[["nsteps"]] <- 365 * 200 
    
    blist[[which(R0vec==R)]] <- replicate(nsim, {
        init <- rinit(pp)
            
        df <- SIRmodel_npatch(pp, init, matrix(1), term_time)
        
        prev <- tail(df$I[df$time%%1==0], 50)/pp[["pop"]]
        
        plot(prev)
        
        data.frame(
            prevalence=prev,
            R0=R
        )
    }, simplify=FALSE)
    
}

bdf <- blist %>%
    lapply(bind_rows) %>%
    bind_rows %>%
    mutate(
        prevalence=prevalence/max(prevalence)
    )

sdf <- reslist %>%
    bind_rows %>%
    group_by(R0, m) %>%
    summarize(
        local=mean(local),
        global=mean(global)
    ) %>%
    gather(key, value, -R0, -m) %>%
    mutate(m=paste0("m = ", m))

gstoch <- ggplot(bdf) +
    geom_point(aes(R0, prevalence), shape=".") +
    geom_line(data=sdf, aes(R0, value, col=key)) +
    scale_y_continuous("Probability", expand=c(0,0)) +
    scale_x_continuous("Reprouctive number", expand=c(0,0)) +
    facet_wrap(~m) +
    theme(
        strip.background = element_blank(),
        panel.spacing = grid::unit(0, "cm"),
        panel.grid = element_blank()
    )

# ggsave("stochastic.pdf", gstoch, width=8, height=6)
