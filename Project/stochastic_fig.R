library(Rcpp)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(directlabels)
if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

L1 <- load("stochastic.rda")
L2 <- load("bifurcation_data.rda")

sdf <- reslist %>%
    bind_rows %>%
    group_by(R0, m) %>%
    summarize(
        local=mean(local),
        global=mean(global)
    ) %>%
    gather(key, value, -R0, -m) %>%
    mutate(m=paste0("m = ", m),
           key=factor(key, levels=c("global", "local"), 
                      labels=c("Global extinction", "Local extinction")))

sdf %>% 
    filter(R0 > 12.5, R0 <=18, m=="m = 0.001") %>%
    group_by(key) %>%
    summarize(
        lwr=range(value)[1],
        upr=range(value)[2]
    )

sdf %>% 
    filter(R0 > 10, key=="Global extinction") %>%
    group_by(m) %>%
    summarize(
        lwr=range(value)[1],
        upr=range(value)[2]
    )
    
bifur_df_norm <- bifur_df %>%
    mutate(lp=log(prevalence)) %>%
    mutate(nlp=lp-min(lp)) %>%
    mutate(nlp=nlp/max(nlp))

gstoch <- ggplot(bifur_df_norm) +
    geom_point(data=filter(bifur_df_norm, i==5, round(R0, 1)==6.8), 
               aes(R0, nlp, group=interaction(i, j), col=factor(i)), size=0.5, alpha=0.5) +
    geom_path(aes(R0, nlp, group=interaction(i, j, sim), col=factor(i)), alpha=0.5) +
    geom_line(data=sdf, aes(R0, value, lty=key), lwd=1.1) +
    scale_y_continuous("Probability of extinction") +
    scale_x_continuous("Reproductive number", expand=c(0,0)) +
    facet_wrap(~m) +
    scale_color_manual(values=c(1, 1, 2, 3, 4, 5), guide=FALSE) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        panel.spacing = grid::unit(0, "cm"),
        panel.grid = element_blank()
    )

ggsave("stochastic.pdf", gstoch, width=7, height=5)

gstoch_ex <- ggplot(bifur_df_norm) +
    geom_point(data=filter(bifur_df_norm, i==5, round(R0, 1)==6.8), 
               aes(R0, nlp, group=interaction(i, j), col=factor(i)), size=0.5, alpha=0.5) +
    geom_path(aes(R0, nlp, group=interaction(i, j, sim), col=factor(i)), alpha=0.5) +
    geom_line(data=filter(sdf, m=="m = 0.001"), aes(R0, value, lty=key), lwd=1.1) +
    geom_dl(data=filter(sdf, m=="m = 0.001"), aes(R0, value, label=key), 
            hjust=2,
            method=list(dl.trans(x= x + .2), "last.bumpup")) +
    scale_y_continuous("Probability of extinction") +
    scale_x_continuous("Reproductive number", expand=c(0,0)) +
    scale_color_manual(values=c(1, 1, 2, 3, 4, 5), guide=FALSE) +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        panel.spacing = grid::unit(0, "cm"),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,6.6,1,1.5), "lines")
    )

gt1 <- ggplotGrob(gstoch_ex)
gt1$layout$clip[gt1$layout$name == "panel"] <- "off"

ggsave("stochastic_example.pdf", gt1, width=7, height=3)
