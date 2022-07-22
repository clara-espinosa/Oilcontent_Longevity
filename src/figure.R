library(tidyverse)

# upload data 
sp <- read.csv("data/species.csv", sep =";")
raw_df <-read.csv("data/germination.csv", sep =";") 
raw_df[is.na(raw_df)] = 0

long_df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  left_join(sp, by="code")%>%
  unite("ecology",distribution:microhabitat, sep=" ", remove = FALSE) %>%
  na.omit  

long_df$success <- long_df$germinated/long_df$seeds 
long_df
View (long_df)

### distribution alpine vs generalist
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

df <- cbind(df, binomials)
 
x11()
ggplot()+
  geom_point(data=long_df, aes(ageing, success, color=distribution), position="jitter", alpha =0.5) +
  scale_color_manual(values=c("dodgerblue4", "orange1")) +
  geom_line(data= df, aes(ageing, mean, color = distribution), size=1.2) +
  geom_errorbar(data= df, aes(ageing, mean, ymin = lower, ymax = upper, color = distribution), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size
  
### Snowbed vs fellfield
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, microhabitat) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

df <- cbind(df, binomials) 
 
ggplot()+
  geom_point(data=long_df, aes(ageing, success, color=microhabitat), position="jitter", alpha =0.5) +
  scale_color_manual(values=c("firebrick2","green3", "deepskyblue")) +
  geom_line(data= df, aes(ageing, mean, color = microhabitat), size=1.2) +
  geom_errorbar(data= df, aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size
        

### Ecolgy alpine/generalist x microhabitatdistribution
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  unite("ecology",distribution:microhabitat, sep=" ", remove = FALSE) %>%
  na.omit %>% 
  group_by(ageing, ecology) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

df <- cbind(df, binomials)
 
ggplot()+
  geom_point(data=long_df, aes(ageing, success, color=ecology), position="jitter", alpha =0.5) +
  scale_color_manual(values=c("firebrick2","green3", "deepskyblue", "yellow")) +
  geom_line(data= df, aes(ageing, mean, color = ecology), size=1.2) +
  geom_errorbar(data= df, aes(ageing, mean, ymin = lower, ymax = upper, color = ecology), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size

### Facet boxplot

raw_df <-read.csv("data/germination.csv", sep =";") 
df <- raw_df %>%
  gather(raw_df, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/species.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, distribution, microhabitat) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df$germinated, df$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat)) +
  facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

cbind(df, binomials) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = distribution)) +
  facet_wrap(~ microhabitat) +
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
