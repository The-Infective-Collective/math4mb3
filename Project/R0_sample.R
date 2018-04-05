library(Rcpp)
library(tidyverse)
source("params.R")
source("useful-fun.R")
sourceCpp("SIRmodel_npatch.cpp")

theme_set(theme_bw())
set.seed(248)

m <- 0.001
M <- matrix(c(1-m, m, m, 1-m), 2, 2)
pp <- base.params
pp[["R0"]] <- 9

init <- r.init.science(2)

df <- run_SIR(pp, init, M, term_time)
#incoherence1 <- apply(df[, c("I1", "I2")], 1, coherence)
  
g1 <- ggplot(df) +
  geom_line(aes(t, I1), color = "#D55E00") +
  geom_line(aes(t, I2), color = "#0072B2") +
  theme_bw() +
  labs(y = "Incidence", x = "Year")
  #geom_line(data = data_frame(t = df$t, incoherence = incoherence1) %, aes(t, incoherence*10000)) +
  #scale_y_continuous(sec.axis = sec_axis(~.*.0001, name = ""))


ggsave("deterministic-sample-coherent.pdf", plot = g1, device = "pdf", width = 12, height = 5, units = "in")

set.seed(2)

init <- r.init.science(2)

df <- run_SIR(pp, init, M, term_time)
#incoherence1 <- apply(df[, c("I1", "I2")], 1, coherence)

g2 <- ggplot(df) +
  geom_line(aes(t, I1), color = "#D55E00") +
  geom_line(aes(t, I2), color = "#0072B2") +
  theme_bw() +
  labs(y = "Incidence", x = "Year")

ggsave("deterministic-sample-incoherent.pdf", plot = g1, device = "pdf", width = 12, height = 5, units = "in")
