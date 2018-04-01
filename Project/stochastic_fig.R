library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

L <- load("stochastic.rda")

sdf <- reslist %>%
    bind_rows %>%
    group_by(R0, m) %>%
    summarize(
        local=mean(local),
        global=mean(global)
    ) %>%
    gather(key, value, -R0, -m) %>%
    mutate(m=paste0("m = ", m))

gstoch <- ggplot(sdf) +
    geom_line(aes(R0, value, col=key)) +
    scale_y_continuous("Probability", expand=c(0,0)) +
    scale_x_continuous("Reprouctive number", expand=c(0,0)) +
    facet_wrap(~m) +
    theme(
        strip.background = element_blank(),
        panel.spacing = grid::unit(0, "cm"),
        panel.grid = element_blank()
    )

# ggsave("stochastic.pdf", gstoch, width=8, height=6)
