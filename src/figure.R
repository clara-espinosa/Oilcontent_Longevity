library(tidyverse)
############## 2023 ######################
# upload data 
sp23 <- read.csv("data/2023/species23.csv", sep =";")
raw23<-read.csv("data/2023/germination23.csv", sep =";") 
raw23[is.na(raw23)] = 0

long23 <- raw23 %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  left_join(sp23, by="code")%>%
  unite("ecology",distribution:microhabitat, sep=" ", remove = FALSE) %>%
  na.omit  

long23$success <- long23$germinated/long23$seeds 
long23
View (long23)

### distribution alpine vs generalist
df23 <- raw23%>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, distribution) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))


binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

df23 <- cbind(df23, binomials23)
 
x11()
ggplot()+
  geom_point(data=long23, aes(ageing, success, color=distribution), position="jitter", alpha =0.5) +
  scale_color_manual(values=c("dodgerblue4", "orange1")) +
  geom_line(data= df23, aes(ageing, mean, color = distribution), size=1.2) +
  geom_errorbar(data= df23, aes(ageing, mean, ymin = lower, ymax = upper, color = distribution), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size
  
### Snowbed vs fellfield
df23 <- raw23%>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, microhabitat) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

df23 <- cbind(df23, binomials23) 
 
ggplot()+
  geom_point(data=long23, aes(ageing, success, color=microhabitat), position="jitter", alpha =0.5) +
  scale_color_manual(values=c("firebrick2","green3", "deepskyblue")) +
  geom_line(data= df23, aes(ageing, mean, color = microhabitat), size=1.2) +
  geom_errorbar(data= df23, aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size
        

### Ecolgy alpine/generalist x microhabitatdistribution
df23 <- raw23 %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  unite("ecology",distribution:microhabitat, sep=" ", remove = FALSE) %>%
  na.omit %>% 
  group_by(ageing, ecology) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

df23 <- cbind(df23, binomials23)
 
ggplot()+
  geom_point(data=long23, aes(ageing, success, color=ecology), position="jitter", alpha =0.5) +
  #scale_color_manual(values=c("firebrick2","green3", "deepskyblue", "yellow")) +
  geom_line(data= df23, aes(ageing, mean, color = ecology), size=1.2) +
  geom_errorbar(data= df23, aes(ageing, mean, ymin = lower, ymax = upper, color = ecology), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size

### Facet boxplot

raw23 <-read.csv("data/2023/germination23.csv", sep =";") 
df23 <- raw23 %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  na.omit %>% 
  group_by(ageing, distribution, microhabitat) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

cbind(df23, binomials23) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat)) +
  facet_wrap(~ distribution) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

cbind(df23, binomials23) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = distribution)) +
  facet_wrap(~ microhabitat) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)

### general boxplots from GerminaR indices
dat23 <- raw23 %>% 
  mutate(across(c(code, ageing), as.factor))
str(dat23)
library (GerminaR)
gerind23 <- ger_summary(SeedN = "seeds", evalName = "D", data=dat23[3:7])
gerind23 <- bind_cols(dat23[,1:2], gerind23)
str(gerind23)

ggplot (gerind23, aes (x=ageing, y=syn, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="syncronization index")+ 
  theme_classic(base_size = 18)

ggplot (gerind23, aes (x=ageing, y=mgr, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="mean germination rate")+ 
  theme_classic(base_size = 18)

ggplot (gerind23, aes (x=ageing, y=grp, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="germination %")+ 
  theme_classic(base_size = 18)


############## 2023 ######################
#based on df23 from mcmc script
raw23 %>%
  mutate(viable = seeds-empty) %>%
  group_by(species, code, ageing) %>%
  mutate(viable = sum(viable)) %>%
  gather(scores, germinated, D0:D62) %>%
  group_by(species,code,  ageing) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE), viable = first(viable)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  rename(seeds = viable) %>%
  mutate (success = germinated/seeds)-> long23


### microhabitat preference
long23 %>%
  group_by(ageing, microhabitat_preference) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds)) ->df23

binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

df23 <- cbind(df23, binomials23)

x11()
ggplot()+
  geom_point(data=long23, aes(ageing, success, color=microhabitat_preference), position="jitter", alpha =0.5) +
  #scale_color_manual(values=c("dodgerblue4", "orange1")) +
  geom_line(data= df23, aes(ageing, mean, color = microhabitat_preference), size=1.2) +
  geom_errorbar(data= df23, aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat_preference), size=1.2) +
  scale_x_continuous (breaks = c(0,2,10,15,30))+
  labs(x = "Ageing days", y="Germination success")+ 
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.8,0.85), 
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="transparent",colour=NA)) #change legend text font size

### Facet boxplot

raw23 <-read.csv("data/2023/longevity23.csv", sep =";") 
df23 <- raw23 %>%
  mutate(viable = seeds-empty) %>%
  group_by(species, code, ageing) %>%
  mutate(viable = sum(viable)) %>%
  gather(scores, germinated, D0:D62) %>%
  group_by(species, code,  ageing) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE), viable = first(viable)) %>%
  merge(read.csv("data/2023/species23.csv", sep =";")) %>%
  rename(seeds = viable) %>%
  group_by(ageing, community, microhabitat_preference) %>%
  summarise(germinated = sum(germinated), seeds = sum(seeds))

binom::binom.confint(df23$germinated, df23$seeds, method = "wilson") %>%
  select(mean:upper) -> binomials23

cbind(df23, binomials23) %>%
  ggplot(aes(ageing, mean, ymin = lower, ymax = upper, color = microhabitat_preference)) +
  facet_wrap(~ community) +
  geom_line(size=1) +
  geom_errorbar(size=1) +
  labs(x = "Ageing days", y="")+ 
  theme_classic(base_size = 18)


### general boxplots from GerminaR indices
dat23 <- raw23 %>% 
  mutate(across(c(code, ageing), as.factor))
str(dat23)
library (GerminaR)
gerind23 <- ger_summary(SeedN = "seeds", evalName = "D", data=dat23[3:7])
gerind23 <- bind_cols(dat23[,1:2], gerind23)
str(gerind23)

ggplot (gerind23, aes (x=ageing, y=syn, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="syncronization index")+ 
  theme_classic(base_size = 18)

ggplot (gerind23, aes (x=ageing, y=mgr, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="mean germination rate")+ 
  theme_classic(base_size = 18)

ggplot (gerind23, aes (x=ageing, y=grp, fill=ageing))+
  geom_boxplot() +
  labs(x = "Ageing days", y="germination %")+ 
  theme_classic(base_size = 18)
