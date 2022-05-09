library(tidyverse)

### Snow bed vs fellfield

read.csv("data/germination.csv") %>%
  gather(days, germinated, D7:D28) %>%
  group_by(code, ageing, sown) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv"), by = "code") %>%
  na.omit %>%
  group_by(ageing, snow) %>%
  summarise(germinated = sum(germinated), sown = sum(sown))-> df

binom::binom.confint(df$germinated, df$sown, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = snow)) +
  # facet_wrap(~ snow) +
  geom_line() +
  geom_errorbar() +
  geom_smooth()
  
### Habitat

read.csv("data/germination.csv") %>%
  gather(days, germinated, D7:D28) %>%
  group_by(code, ageing, sown) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv"), by = "code") %>%
  na.omit %>%
  group_by(ageing, habitat) %>%
  summarise(germinated = sum(germinated), sown = sum(sown))-> df

binom::binom.confint(df$germinated, df$sown, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = habitat)) +
  # facet_wrap(~ habitat) +
  geom_line() +
  geom_errorbar() +
  geom_smooth()

### Bedrock

read.csv("data/germination.csv") %>%
  gather(days, germinated, D7:D28) %>%
  group_by(code, ageing, sown) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv"), by = "code") %>%
  na.omit %>%
  group_by(ageing, bedrock) %>%
  summarise(germinated = sum(germinated), sown = sum(sown))-> df

binom::binom.confint(df$germinated, df$sown, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = bedrock)) +
  # facet_wrap(~ bedrock) +
  geom_line() +
  geom_errorbar() +
  geom_smooth()
