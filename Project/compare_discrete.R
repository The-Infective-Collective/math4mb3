library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(deSolve)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

SIR <- function(t, y, param) {
    with(as.list(c(y, param)), {
        beta <- term_time(t, param)
        
        dS <- mu * (pop-S) - beta * S * I
        dI <- beta * S * I - (mu + gamma) * I
        dR <- gamma * I - mu * R
        
        list(c(dS, dI, dR))
    })
}

nyear <- 40

base.params$nsteps <- 365 * nyear

tvec <- seq(0, nyear, by=1/365)

reslist <- vector('list', 3)

ii <-list(
    S=10000,
    I=2000,
    R=1e6-10000-2000
)


pp <- base.params
    
cont.df <- as.data.frame(deSolve::ode(unlist(ii), tvec, SIR, base.params))
disc.df <- SIRmodel_npatch(base.params, ii, matrix(1), term_time)
    
sl <- list(
    continuous=cont.df[,c("time", "I")],
    discrete=data.frame(
        time=disc.df$time,
        I=disc.df$I
    )
) %>%
    bind_rows(.id="type")

rdf <- sl %>%
    bind_rows(.id="sim") %>%
    mutate(
        sim=factor(sim, levels=c(1,2), 
                   labels=c("S(0)=5000, I(0)=1000",
                            "S(0)=10000, I(0)=2000"))
    )

gcompare <- ggplot(rdf) +
    geom_line(aes(time, I, col=type, linetype=type),lwd=1.1) +
    scale_x_continuous("Time (years)") +
    scale_y_continuous("Prevalence") +
    scale_color_manual(values=c("#D55E00", "#0072B2")) +
    scale_linetype_manual(values=c(1,2)) +
    facet_wrap(~sim, scale="free", ncol=1) +
    theme(
        strip.background = element_blank(),
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.9),
        legend.key.width = grid::unit(2, "cm")
    )

ggsave("compare.pdf", gcompare, width=8, height=5)
