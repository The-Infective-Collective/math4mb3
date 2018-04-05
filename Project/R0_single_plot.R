library(tidyverse)
library(latex2exp)

# load data correctly
load("R0m.rda")
load("bifurcation_data.rda")
df1 <- reslist

R0df <- bind_rows(df1[[1]]) %>%
  group_by(R0, m) %>%
  summarize(prob.coherence = length(which((incoherence2 < 100)))/100)


bifur_df_norm <- bifur_df %>%
  mutate(lp=log(prevalence)) %>%
  mutate(nlp=lp-min(lp)) %>%
  mutate(nlp=nlp/max(nlp))


gprob <- ggplot(R0df) +
  geom_line(aes(R0, prob.coherence), size=1) +
  geom_path(data= bifur_df_norm, aes(R0, nlp, group=interaction(i, j, sim), col=factor(i)), alpha=0.5) +
  labs(x='Reproductive number', y='Probability of coherence') + 
  theme_bw()+
  scale_color_manual(values=c(1,1,2,3,4,5,6))

ggsave("probabilitycoherence10-3.pdf", gprob, width=12, height=5)
