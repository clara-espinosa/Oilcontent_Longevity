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
  #facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)
 # geom_smooth()
 

### distribution alpine vs generalist
raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, distribution) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = distribution)) +
  # facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)
  

### Ecolgy alpine/generalist x microdistribution
raw_df <- read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  unite("ecology",distribution:micro, sep=" ", remove = FALSE) %>%
  na.omit %>% 
  group_by(ageing, ecology) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = ecology)) +
  # facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

### Facet boxplot

raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, distribution, micro) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = micro)) +
  facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = distribution)) +
  facet_wrap(~ micro) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

### general boxplots from GerminaR indices
dat <- raw_df %>% 
  mutate(across(c(code, ageing), as.factor))
str(dat)
library (GerminaR)
gerind <- ger_summary(SeedN = "seeds", evalName = "D", data=dat[3:7])
gerind <- bind_cols(dat[,1:2], gerind)
str(gerind)

ggplot (gerind, aes (x=ageing, y=syn, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="syncronization index")+ 
  theme_classic(base_size = 18)

ggplot (gerind, aes (x=ageing, y=mgr, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="mean germination rate")+ 
  theme_classic(base_size = 18)

ggplot (gerind, aes (x=ageing, y=grp, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="germination %")+ 
  theme_classic(base_size = 18)
