library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

L <- load("bifurcation.rda")

save <- FALSE

bdf <- blist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows

bifur_df <- bdf %>%
    group_by(sim, R0) %>%
    filter(prevalence > 0) %>%
    filter(!duplicated(prevalence)) %>%
    group_by(R0) %>%
    filter(!duplicated(round(prevalence, 7))) %>%
    group_by(sim, R0) %>%
    mutate(
        i=n()
    ) %>%
    group_by(R0, i) %>%
    mutate(
        prevalence=ifelse(i==1, mean(prevalence), prevalence)
    ) %>%
    filter(!duplicated(prevalence)) %>%
    group_by(sim, R0, i) %>%
    mutate(
        prevalence=sort(prevalence),
        j=seq_along(i)
    ) %>%
    group_by(R0, i, j) %>%
    summarize(
        prevalence=mean(prevalence)
    ) %>%
    filter(
        !(i==2 & R0 < 13.6),
        !(i==1 & R0 > 13.8),
        !(i==3 & R0 > 15.5 & R0 < 15.7)
    ) %>%
    as.data.frame

gbifur <- ggplot(bifur_df) +
    geom_path(aes(R0, prevalence, group=interaction(i, j), col=factor(i)), lwd=2) +
    scale_y_log10("Prevalence I/N", breaks=c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)) +
    scale_x_continuous("Basic reproductive number") +
    scale_color_manual(values=c(1, 1, 2, 3, 4, 5)) +
    theme(
        legend.position = "none",
        panel.grid = element_blank()
    )

if (save) ggsave("bifurcation.pdf", gbifur, width=10, height=6)

if (save) save("bifur_df", file="bifurcation_data.rda")
