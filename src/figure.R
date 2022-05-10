library(tidyverse)

### Snowbed vs fellfield
raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, micro) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = micro)) +
  #facet_wrap(~ habitat) +
  geom_line() +
  geom_errorbar() 
 # geom_smooth()
 
### Habitat
raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, habitat) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = habitat)) +
  # facet_wrap(~ habitat) +
  geom_line() +
  geom_errorbar() 
  # geom_smooth()


