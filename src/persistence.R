library(tidyverse);library (rstatix);library (stringr)
library(ggpattern); library (vegan);library(binom)
library(glmmTMB); library (DHARMa)

##### PERSISTENCE 23 #####
### Read data
library(viridis)
str(persistence)
x11()
# work in progress not clear if field germination should be
# subset of germination responses on the lab that will be used to compare both microclimatic
# conditions
persistence%>%
  select(species, community, retrieval_season, site_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  filter(retrieval_season=="Summer_22")%>%
  mutate (micro1 = "Warm", 
          micro2="Cold")%>%
  gather(micro,microhabitat_buried, micro1:micro2)%>%
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)-> lab_duplicate
  
  
# data handling and ggplot with point and error bars x species
persistence %>%
  filter(!retrieval_season=="Summer_22")%>% # remove original lab data to then add the duplicated
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  rbind(lab_duplicate)%>%  
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  filter (!retrieval_season == "Spring_24")%>%
  filter (!retrieval_season == "Autumn_24")%>%
  replace(is.na(.), 0)%>%
  rename(Field_germ=D0)%>%
  mutate(labgerm = rowSums(across(D7:D98)))%>%
  select(species, community, site_buried, microhabitat_buried, retrieval_season, 
         seeds_initial, seeds_identifiables, Field_germ, labgerm, empty, fungi)%>%
  group_by(species, community, site_buried, microhabitat_buried, retrieval_season)%>%
  mutate(not_persistent = seeds_initial-seeds_identifiables)%>%
  mutate (no_viable = fungi+empty) %>%
  mutate(nogermbutviable = seeds_identifiables - (labgerm+Field_germ+no_viable))%>%
  mutate (remain_viable = labgerm+nogermbutviable) %>%
  group_by (community, species, microhabitat_buried, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_identifiables), 
            seeds_initial = sum (seeds_initial))%>%
  mutate (binom.confint(remain_viable, seeds_initial, methods = "wilson")) %>% 
  filter(community =="mediterranean") %>% #mediterranean, temperate
  mutate(retrieval_season = as.factor(retrieval_season))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Summer_22","Spring_23","Autumn_23"))%>%
  ggplot(aes(x=retrieval_season, y=mean, fill= microhabitat_buried)) +
  geom_point(shape=21, size =4) + #, position = "jitter"
  geom_errorbar(aes(retrieval_season, mean, ymin = lower, ymax = upper, color = microhabitat_buried),width = 0.2, size =1) +
  scale_fill_manual (name= "Microhabitat buried",values = c ("chocolate2","deepskyblue3")) + #, guide = "none"
  geom_line(aes(x=retrieval_season, y=mean, group= microhabitat_buried, color = microhabitat_buried), linewidth = 1)+
  scale_color_manual (values = c ("chocolate2","deepskyblue3"), guide = "none") +
  facet_wrap(~species, ncol= 2)+
  labs (title = "Mediterranean species", x= "Retrieval season", y= "Mean Seed Viability")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text (size = 18),
        strip.text= element_text(face="italic", size=14),
        legend.position= "bottom")

# data handling and ggplot with point and error bars x all mediterranean
persistence %>%
  filter(!retrieval_season=="Summer_22")%>% # remove original lab data to then add the duplicated
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  rbind(lab_duplicate)%>%  
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  filter (!retrieval_season == "Spring_24")%>%
  filter (!retrieval_season == "Autumn_24")%>%
  replace(is.na(.), 0)%>%
  rename(Field_germ=D0)%>%
  mutate(labgerm = rowSums(across(D7:D98)))%>%
  select(species, community, site_buried, microhabitat_buried, retrieval_season, 
         seeds_initial, seeds_identifiables, Field_germ, labgerm, empty, fungi)%>%
  group_by(species, community, site_buried, microhabitat_buried, retrieval_season)%>%
  mutate(not_persistent = seeds_initial-seeds_identifiables)%>%
  mutate (no_viable = fungi+empty) %>%
  mutate(nogermbutviable = seeds_identifiables - (labgerm+Field_germ+no_viable))%>%
  mutate (remain_viable = labgerm+nogermbutviable) %>%
  group_by (community, species, microhabitat_buried, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_identifiables), 
            seeds_initial = sum (seeds_initial))%>%
  group_by (community, microhabitat_buried, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_id), 
            seeds_initial = sum (seeds_initial))%>%
  mutate (binom.confint(remain_viable, seeds_initial, methods = "wilson")) %>% 
  #filter(community =="mediterranean") %>% #mediterranean, temperate
  mutate(retrieval_season = as.factor(retrieval_season))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Summer_22","Spring_23","Autumn_23"))%>%
  ggplot(aes(x=retrieval_season, y=mean, fill= microhabitat_buried)) +
  geom_point(shape=21, size =4) + #, position = "jitter"
  geom_errorbar(aes(retrieval_season, mean, ymin = lower, ymax = upper, color = microhabitat_buried),width = 0.2, size =1) +
  scale_fill_manual (name = "Microhabitat buried",values = c ("chocolate2","deepskyblue3")) + #, guide = "none"
  geom_line(aes(x=retrieval_season, y=mean, group= microhabitat_buried, color = microhabitat_buried), linewidth = 1)+
  scale_color_manual (values = c ("chocolate2","deepskyblue3"), guide = "none") +
  facet_wrap(~community, ncol= 2)+
  #scale_color_viridis_d() +
  #ylim (c(0,1))+
  labs (title = "Viability Loss per Community", x= "Retrieval season", y= "Mean Seed Viability")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text (size = 18),
        strip.text= element_text(face="bold", size=14),
        legend.position= "bottom")

# data handling and ggplot with stacked bar (remain viable, field germ and (no germ+ no persistent)
library(RColorBrewer)
x11()
persistence %>%
  filter(!retrieval_season=="Summer_22")%>% # remove original lab data to then add the duplicated
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  rbind(lab_duplicate)%>%  
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  filter (!retrieval_season == "Spring_24")%>%
  filter (!retrieval_season == "Autumn_24")%>%
  replace(is.na(.), 0)%>%
  rename(Field_germ=D0)%>%
  mutate(labgerm = rowSums(across(D7:D98)))%>%
  select(species, community, site_buried, microhabitat_buried, retrieval_season, 
         seeds_initial, seeds_identifiables, Field_germ, labgerm, empty, fungi)%>%
  group_by(species, community, site_buried, microhabitat_buried, retrieval_season)%>%
  mutate(not_persistent = seeds_initial-seeds_identifiables)%>%
  mutate (no_viable = fungi+empty) %>%
  mutate(nogermbutviable = seeds_identifiables - (labgerm+Field_germ+no_viable))%>%
  mutate (remain_viable = labgerm+nogermbutviable) %>%
  group_by (community, species, microhabitat_buried, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_identifiables), 
            seeds_initial = sum (seeds_initial))%>%
  mutate(no_viable = no_viable+not_persistent)%>%
  select(community, species, microhabitat_buried, retrieval_season, Field_germ, remain_viable, no_viable)%>%
  gather(seed_status, num_seeds,Field_germ:no_viable)%>%
  #group_by (community, microhabitat_buried, retrieval_season)%>%
  #summarise(Field_germ = sum(Field_germ), 
  #          remain_viable = sum(remain_viable), 
  #          no_viable= sum(no_viable),
  #          not_persistent=sum(not_persistent),
  #          seeds_id = sum(seeds_id), 
  #          seeds_initial = sum (seeds_initial))%>%
  #mutate (binom.confint(remain_viable, seeds_initial, methods = "wilson")) %>% 
  filter(community =="temperate") %>% #mediterranean, temperate
  mutate(retrieval_season = as.factor(retrieval_season))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Summer_22","Spring_23","Autumn_23"))%>%
  ggplot(aes(x=microhabitat_buried, y=num_seeds, fill= seed_status)) +
  geom_bar(position="fill", stat="identity") + #, position = "stack" = count; "fill" trnasform to percentages
  #geom_errorbar(aes(retrieval_season, mean, ymin = lower, ymax = upper, color = microhabitat_buried),width = 0.2, size =1) +
  scale_fill_manual (name = "Seed status",values = c ("gold","chocolate","darkred")) + #, guide = "none"
  #geom_line(aes(x=retrieval_season, y=mean, group= microhabitat_buried, color = microhabitat_buried), linewidth = 1)+
  #scale_color_manual (values = c ("deepskyblue3","chocolate2"), guide = "none") +
  facet_wrap(~retrieval_season)+#, ncol= 2
  #scale_color_viridis_d() +
  #ylim (c(0,1))+
  labs (title = "Temperate", x= "Retrieval season", y= "Seeds proportion")+
  theme_classic(base_size = 12)+
  theme(plot.title = element_text (size = 18),
        strip.text= element_text(face="italic", size=14),
        legend.position= "bottom")


##### Field germination for move-along graph + glms ####

x11()
#bar plot
field_germ %>%
  filter(retrieval_season == "Spring_23" |  retrieval_season == "Autumn_23")%>%
  select(retrieval_season, community, microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, microhabitat_buried, species) %>%
  mutate(species = fct_relevel (species, "Armeria duriaei", "Dianthus langeanus", "Plantago holosteum",
                                "Luzula caespitosa", "Phyteuma hemisphaericum", "Silene ciliata", 
                                "Androsace villosa",  "Carex sempervirens","Gypsophila repens", 
                                "Armeria cantabrica","Festuca glacialis","Jasione cavanillesii"))%>%
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  gather(bag, germ, bag1:bag3)%>%
  #mutate(germ_per=germ/10)%>%
  group_by(retrieval_season, microhabitat_buried, species)%>%
  summarise(germ=sum(germ))%>%
  #print(n=48)
  spread(retrieval_season, germ)%>%
  mutate(Autumn_23 = Autumn_23-Spring_23)%>%
  mutate(Autumn_23 = ifelse(Autumn_23>0,Autumn_23, 0)) %>%
  gather(retrieval_season, field_germ, Spring_23:Autumn_23)%>%
  mutate(seeds_initial = 60)%>%
  mutate (binom.confint(field_germ, seeds_initial, methods = "wilson")) %>% 
  mutate(retrieval_season = recode_factor (retrieval_season, "Spring_23" = "Spring","Autumn_23"= "Autumn"))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Spring","Autumn"))%>%
  #filter(species=="Armeria duriaei")%>%
  ggplot(aes(retrieval_season, mean, condition=microhabitat_buried, fill = microhabitat_buried))+
  geom_bar(width = 0.7, position=position_dodge(width = 0.8), stat="identity", color="black")+
  geom_errorbar(aes(retrieval_season, mean, ymin = lower, ymax = upper), color="black",position=position_dodge(width = 0.8) ,width = 0.2, size =1, show.legend = F) +
  scale_fill_manual (name="Microhabitat",values = c ("chocolate2", "deepskyblue3")) +
  scale_color_manual (values = c ("chocolate2", "deepskyblue3")) + #
  labs (y= "Field germinated seeds", x= "Retrieval season")+ #title= "Field germination",
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  facet_wrap(~species, ncol=3)+
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        strip.text.x = element_text( size = 12, face = "italic"),# face = "bold",
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "bottom", #bottom
        plot.tag.position = c(0,1),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=15), 
        axis.title.x = element_text (size=15))

#dataframe for testing field differences
field_germ %>%
  filter(retrieval_season == "Spring_23" |  retrieval_season == "Autumn_23")%>%
  select(retrieval_season, community, site_buried,microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, site_buried, microhabitat_buried, species) %>%
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  gather(bag, germ, bag1:bag3)%>%
  mutate(seeds_initial = 10) %>%
  #filter(retrieval_season == "Spring_23" ) -> spring_df
  filter(retrieval_season == "Autumn_23")-> autumn_df



### Get GLM coefficients

glms <- function(x) {
  glm(cbind(germ, seeds_initial - germ) ~ microhabitat_buried, family = "binomial", data = x) -> m1 # site_buried only 2 levels, 
  broom::tidy(m1)
}

autumn_df%>%
  group_by(species) %>%
  do(glms(.))%>%
  print(n=36)
#### IBC Persistence ####
# data handling and ggplot with stacked bar (remain viable, field germ and (no germ+ no persistent)
library(RColorBrewer)
x11()
persistence %>%
  filter(!retrieval_season=="Summer_22")%>% # remove original lab data to then add the duplicated
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  rbind(lab_duplicate)%>%  
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  filter (!retrieval_season == "Spring_24")%>%
  filter (!retrieval_season == "Autumn_24")%>%
  replace(is.na(.), 0)%>%
  rename(Field_germ=D0)%>%
  mutate(labgerm = rowSums(across(D7:D98)))%>%
  select(species, community, site_buried, microhabitat_buried, retrieval_season, 
         seeds_initial, seeds_identifiables, Field_germ, labgerm, empty, fungi)%>%
  group_by(species, community, site_buried, microhabitat_buried, retrieval_season)%>%
  mutate(not_persistent = seeds_initial-seeds_identifiables)%>%
  mutate (no_viable = fungi+empty) %>%
  mutate(nogermbutviable = seeds_identifiables - (labgerm+Field_germ+no_viable))%>%
  mutate (remain_viable = labgerm+nogermbutviable) %>%
  group_by (community, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_identifiables), 
            seeds_initial = sum (seeds_initial))%>%
  mutate(no_viable = no_viable+not_persistent)%>%
  select(community, retrieval_season, Field_germ, remain_viable, no_viable)%>%
  gather(seed_status, num_seeds,Field_germ:no_viable)%>%
  mutate(retrieval_season = as.factor(retrieval_season))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Summer_22","Spring_23","Autumn_23"))%>%
  mutate(community = as.factor(community))%>%
  mutate(community = recode (community, "temperate"="Temperate" , "mediterranean"= "Mediterranean"))%>%
  #mutate(community= fct_relevel(community, "Temperate","Mediterranean")) %>%
  ggplot(aes(x=retrieval_season, y=num_seeds, fill= seed_status), color="black") +
  geom_bar(position="fill", stat="identity", color="black") + #, position = "stack" = count; "fill" trnasform to percentages
  scale_fill_manual (name = "Seed status",values = c ("chocolate2","darkgrey","deepskyblue3")) + #, guide = "none"
  facet_wrap(~community)+#, ncol= 2
  geom_hline(yintercept=0.5, linetype ="dashed", size =1.5, colour = "black") +
  labs (x= "Retrieval season", y= "Seeds proportion")+
  theme_classic(base_size = 16)+
  theme(plot.title = element_text (size = 18),
        strip.text= element_text(size=20),
        legend.position= "bottom")

# IBC GLM
persistence %>%
  filter(!retrieval_season=="Summer_22")%>% # remove original lab data to then add the duplicated
  select(species, community, retrieval_season, site_buried, microhabitat_buried, 
         seeds_initial, seeds_identifiables, D0:fungi)%>%
  rbind(lab_duplicate)%>%  
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  filter (!retrieval_season == "Spring_24")%>%
  filter (!retrieval_season == "Autumn_24")%>%
  replace(is.na(.), 0)%>%
  rename(Field_germ=D0)%>%
  mutate(labgerm = rowSums(across(D7:D98)))%>%
  select(species, community, site_buried, microhabitat_buried, retrieval_season, 
         seeds_initial, seeds_identifiables, Field_germ, labgerm, empty, fungi)%>%
  group_by(species, community, site_buried, microhabitat_buried, retrieval_season)%>%
  mutate(not_persistent = seeds_initial-seeds_identifiables)%>%
  mutate (no_viable = fungi+empty) %>%
  mutate(nogermbutviable = seeds_identifiables - (labgerm+Field_germ+no_viable))%>%
  mutate (remain_viable = labgerm+nogermbutviable) %>%
  group_by (community, species, retrieval_season)%>%
  summarise(Field_germ = sum(Field_germ), 
            remain_viable = sum(remain_viable), 
            no_viable= sum(no_viable),
            not_persistent=sum(not_persistent),
            seeds_id = sum(seeds_identifiables), 
            seeds_initial = sum (seeds_initial))%>%
  mutate(no_viable = no_viable+not_persistent)%>%
  select(community, species, retrieval_season, Field_germ, remain_viable, no_viable,  seeds_initial)%>%
  #merge(read.csv("data/2023/species23.csv", sep =",")) %>%
  #select(community, species, retrieval_season, Field_germ, remain_viable, no_viable)%>%
  mutate(ID = gsub(" ", "_", species), animal = ID)-> IBC
  #select(!(family))

### GLM without phylogeny

glm(cbind(remain_viable, seeds_initial - remain_viable) ~
      community*retrieval_season, family = "binomial", data = IBC) -> m2
summary(m2)

#EXTRA VISUALIZATION #####
#faceted box plot
field_germ %>%
  filter(retrieval_season == "Spring_23" |  retrieval_season == "Autumn_23")%>%
  select(retrieval_season, community, microhabitat_buried, species, seeds_initial, bag1:bag3)%>%
  convert_as_factor(retrieval_season, community, microhabitat_buried, species) %>%
  mutate(species = fct_relevel (species, "Armeria duriaei", "Dianthus langeanus", "Plantago holosteum",
                                "Luzula caespitosa", "Phyteuma hemisphaericum", "Silene ciliata", 
                                "Androsace villosa",  "Carex sempervirens","Gypsophila repens", 
                                "Armeria cantabrica","Festuca glacialis","Jasione cavanillesii"))%>%
  mutate(retrieval_season = fct_relevel (retrieval_season, "Spring_23","Autumn_23"))%>%
  mutate(microhabitat_buried = recode_factor(microhabitat_buried, "Warm" = "Fellfield", 
                                             "Crio"="Fellfield", 
                                             "Cold" = "Snowbed", 
                                             "Snow" = "Snowbed"))%>%
  gather(bag, germ, bag1:bag3)%>%
  mutate(germ_per=germ/10)%>%
  filter(species=="Armeria duriaei")%>%
  ggplot(aes(microhabitat_buried, germ_per, pattern=retrieval_season, fill=microhabitat_buried))+
  geom_boxplot(color="black")+
  geom_boxplot_pattern(position = position_dodge(preserve = "single"), color = "black", 
                       pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, 
                       pattern_spacing = 0.05, pattern_key_scale_factor = 0.6) +
  scale_fill_manual (values = c ("chocolate2", "deepskyblue3"), guide = "none") + #
  labs (y= "Germination proportion", x= "Microhabitat buried")+ #title= "Field germination",
  scale_y_continuous (limits = c(0,1), breaks = seq (0, 1, by= 0.25)) +
  facet_wrap(~species, ncol=3)+
  theme_classic(base_size = 14) +
  theme(plot.title = element_text (size = 22),
        strip.text.x = element_text( size = 14, face = "italic"),# face = "bold",
        strip.text.y = element_text(size = 14, angle = 360),
        legend.position = "right", #bottom
        plot.tag.position = c(0,1),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text (size=15), 
        axis.title.x = element_text (size=15))