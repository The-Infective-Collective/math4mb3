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

## LOCAL EXTINCTION

pp <- base.params
pp[["R0"]] <- 8.4

init <- initfun(pp, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df2 <- SIRmodel_npatch_stochastic(pp, init, M, term_time)

simdf1 <- list(
    patch1=data.frame(
        time=df2$time,
        I=df2$I[,1]
    ),
    patch2=data.frame(
        time=df2$time,
        I=df2$I[,2]
    )
) %>%
    bind_rows(.id="patch") %>%
    gather(key, value, -time, -patch) %>%
    filter(time > 20, time < 40)

extdf1 <-simdf1 %>%
    filter(value==0)

gbase <- ggplot(simdf1) +
    geom_line(aes(time, value, col=patch), lwd=1.1) +
    labs(x="Time (years)", y="Prevalence")+
    ggtitle("Local extinction (with rescue effects)") +
    scale_color_manual(labels=c("Patch 1", "Patch 2"), values=c("#D55E00", "#0072B2")) +
    theme(
        legend.position = c(0.15, 0.9),
        legend.title = element_blank()
    )

glocal <- gbase + geom_point(data=extdf1, aes(time, value), col="red")

# ggsave("localext_1.pdf", glocal, width=8, height=6)

## GLOBAL EXTINCTION

pp2 <- base.params
pp2[["R0"]] <- 2.0

init <- initfun(pp2, 2, T)
m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df3 <- SIRmodel_npatch_stochastic(pp2, init, M, term_time)

simdf2 <- list(
    patch1=data.frame(
        time=df3$time,
        I=df3$I[,1]
    ),
    patch2=data.frame(
        time=df3$time,
        I=df3$I[,2]
    )
) %>%
    bind_rows(.id="patch") %>%
    gather(key, value, -time, -patch) %>%
    filter(time > 20, time < 40)

#extdf2 <-simdf2 %>%
    #filter(value==0)

gglobal <- (gbase %+% simdf2) +
    ggtitle("Global extinction (no rescue effects)") +
    theme(
        legend.position = "none"
    )
    
# ggsave("globalext_1.pdf", gglobal, width=8, height=6)

gext <- arrangeGrob(glocal, gglobal, nrow=1)

ggsave("extinction.pdf", gext, width=10, height=6)

## ASYNCHRONY WITH LOW M

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
    filter(time > 0, time < 20)

extdf3 <-simdf3 %>%
    filter(value==0)

gbase2 <- ggplot(simdf3) +
    geom_line(aes(time, value, col=patch)) +
    labs(x="Time (years)", y="Prevalence")+
    ggtitle("Observed asynchrony with low mixing rate (m = 0.001)")+
    scale_color_manual(labels=c("Patch 1", "Patch 2"), values=c(1,2)) +
    theme(
        legend.position = c(0.05, 0.8),
        legend.title = element_blank()
    )

gasynch <- gbase2 + geom_point(data=extdf3, aes(time, value), col="blue")

#ggsave("asynchrony_lowm_1.pdf", gasynch, width=8, height=6)

## SYNCHRONY WITH HIGH M

pp4 <- base.params
pp4[["R0"]] <- 17.0

init <- initfun(pp4, 2, T)
m <- 0.5
M <- matrix(c(1-m, m, m, 1-m), 2, 2)

set.seed(10)
df5 <- SIRmodel_npatch_stochastic(pp4, init, M, term_time)

simdf4 <- list(
    patch1=data.frame(
        time=df5$time,
        I=df5$I[,1]
    ),
    patch2=data.frame(
        time=df5$time,
        I=df5$I[,2]
    )
) %>%
    bind_rows(.id="patch") %>%
    gather(key, value, -time, -patch) %>%
    filter(time > 0, time < 20)

extdf4 <-simdf4 %>%
    filter(value==0)

gbase3 <- ggplot(simdf4) +
    geom_line(aes(time, value, col=patch)) +
    labs(x="Time (years)", y="Prevalence")+
    ggtitle("Observed synchrony with high mixing rate (m = 0.5)")+
    scale_color_manual(labels=c("Patch 1", "Patch 2"), values=c(1,2)) +
    theme(
        legend.position = "none"
        )

gsynch <- gbase3 + geom_point(data=extdf4, aes(time, value), col="blue")

gtraj <- arrangeGrob(gasynch, gsynch, ncol=1)

ggsave("trajectories.pdf", gtraj, width=10, height=6)
