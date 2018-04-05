library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

beta_df <- data.frame(
    t=t <- seq(0, 1, by=1/3650),
    beta=sapply(t, term_time, base.params)
)

ggterm <- ggplot(beta_df) +
    geom_line(aes(t, beta)) +
    scale_y_continuous(expression(beta(t))) +
    scale_x_continuous("Time (years)") +
    theme(
        panel.grid = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    )
ggsave("term_time.pdf", ggterm, width=8, height=3)




