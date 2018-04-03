library(tidyverse)

# load data correctly
load("R0m.rda")
df1 <- reslist

load("R0.rda")
df2 <- reslist

R0list <- c(df1, df2)
R0df <- do.call("rbind", R0list)

df <- R0df %>%
  group_by(R0, m) %>%
  summarize(prob.coherence = length(which((incoherence2 < 0.2)))/100)

ggplot(df) +
  geom_line(aes(R0, prob.coherence)) +
  facet_wrap(~m, ncol =2)
