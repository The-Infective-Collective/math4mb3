library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

library(Rcpp)
library(gridExtra)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

pp3 <- base.params
pp3[["R0"]] <- 17.0

init <- initfun(pp3, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df4 <- SIRmodel_npatch_stochastic(pp3, init, M, term_time)

simdf3 <- list(
    patch1=data.frame(
        time=df4$time,
        I=df4$I[,1]
    ),
    patch2=data.frame(
        time=df4$time,
        I=df4$I[,2]
    )
) %>%
    bind_rows(.id="patch") %>%
    gather(key, value, -time, -patch) %>%
    filter(time > 60, time < 85)

gex <- ggplot(simdf3) +
    geom_line(aes(time, value, col=patch), lwd=1.1) +
    scale_x_continuous("Time (years)", expand=c(0, 0)) +
    scale_y_continuous("Prevalence") +
    scale_color_manual(labels=c("Patch 1", "Patch 2"), values=c(1, 2)) +
    theme(
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.05, 0.85)
    )

ggsave("stochastic_illustrate.pdf", gex, width=10, height=4)

