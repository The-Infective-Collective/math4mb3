library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(Rcpp)
source("params.R")
sourceCpp("SIRmodel_npatch.cpp")

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

glocal <- ggplot(simdf1) +
    geom_line(aes(time, value, col=patch)) +
    geom_point(data=extdf1, aes(time, value), col=3)+
    labs(x="Time (years)", y="Prevalence")

ggsave("localext_1.pdf", glocal, width=8, height=6)

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

extdf2 <-simdf2 %>%
    filter(value==0)

gglobal <- ggplot(simdf2) +
    geom_line(aes(time, value, col=patch)) +
    geom_point(data=extdf2, aes(time, value), col=3)+
    labs(x="Time (years)", y="Prevalence")

ggsave("globalext2.pdf", gglobal, width=8, height=6)

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
    filter(time > 0, time < 50)

extdf3 <-simdf3 %>%
    filter(value==0)

gasynch <- ggplot(simdf3) +
    geom_line(aes(time, value, col=patch)) +
    geom_point(data=extdf3, aes(time, value), col=3)+
    labs(x="Time (years)", y="Prevalence")

ggsave("asynchrony_lowm2.pdf", gasynch, width=8, height=6)

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
    filter(time > 0, time < 50)

extdf4 <-simdf4 %>%
    filter(value==0)

gsynch <- ggplot(simdf4) +
    geom_line(aes(time, value, col=patch)) +
    geom_point(data=extdf4, aes(time, value), col=3)+
    labs(x="Time (years)", y="Prevalence")

ggsave("synchrony_highm2.pdf", gsynch, width=8, height=6)
