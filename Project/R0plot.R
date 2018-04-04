library(tidyverse)
library(latex2exp)

# load data correctly
load("R0m.rda")
load("bifurcation_data.rda")
df1 <- reslist

load("R0.rda")
df2 <- reslist

R0df <- bind_rows(df1,df2) %>%
    group_by(R0, m) %>%
    summarize(prob.coherence = length(which((incoherence2 < 0.2)))/100)


bifur_df_norm <- bifur_df %>%
    mutate(lp=log(prevalence)) %>%
    mutate(nlp=lp-min(lp)) %>%
    mutate(nlp=nlp/max(nlp))


ggplot(R0df) +
  geom_line(aes(R0, prob.coherence)) +
  facet_wrap(~m, ncol =2) + 
  geom_path(data= bifur_df_norm, aes(R0, nlp, group=interaction(i, j, sim), col=factor(i)), alpha=0.5) +
  labs(x=TeX('$R_0$'), y='Probability of coherence')

