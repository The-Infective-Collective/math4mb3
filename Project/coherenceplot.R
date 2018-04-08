library(tidyverse)
library(latex2exp)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                   base_family = "Times"))

load("bifurcation_data.rda")

load("coherence_m0.0001.rda")
df1 <- reslist

load("coherence_m0.001-0.01.rda")
df2 <- reslist

load("coherence_m0.1-0.5.rda")
df3 <- reslist


R0df <- bind_rows(df1, df2,df3) %>% 
    group_by(R0, m) %>%
    summarize(prob.coherence = length(which((incoherence2 < 100)))/100) %>%
    mutate(m.factor = as.factor(m))

bifur_df_norm <- bifur_df %>%
    mutate(lp=log(prevalence)) %>%
    mutate(nlp=lp-min(lp)) %>%
    mutate(nlp=nlp/max(nlp))

levels(R0df$m.factor) <- c("m = 0.0001", "m = 0.001", "m = 0.01", "m = 0.1", "m = 0.5")

gprob <- ggplot(R0df) +
  geom_line(aes(R0, prob.coherence),size=1) +
  facet_wrap(~m.factor, ncol =2) + 
  geom_path(data= bifur_df_norm, aes(R0, nlp, group=interaction(i, j, sim), col=factor(i)), alpha=0.5) +
  labs(x='Reproductive number', y='Probability of coherence') + 
  scale_color_manual(values=c(1,1,2,3,4,5,6))+
  theme(legend.position="none")

ggsave("probabilitycoherence.png", device="png", gprob, width=7, height=5)