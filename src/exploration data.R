library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan)

# oil data composition and clasification according to n carbons
read.csv("data/oil_type_specification.csv")

# species explanatory data
read.csv("data/species.csv")%>%
  get_summary_stats()
# exploration graph of oil content data

#### 1- Total oil content (%) per seed ####
# 1.1 oil content per sp (mg oil /g sample)
read.csv("data/oil_data_long.csv")%>%
  group_by(Taxon, species, family) %>%
  summarise(PERoil = mean(PERoil))%>%
  ggplot(aes(x=Taxon, y= PERoil, fill=species ))+ #community
  #geom_point(size = 4, shape=21, color = "black", show.legend = F)+
  geom_bar(stat = "identity", show.legend = F )+ #problem with Thymus and S. ciliata (2 accesions summ % oil)
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "Species oil content (%)", tag = "", y= "Oil content (%)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10, color = "black"))
  
# 1.2 oil content per family (mg oil /g sample)
read.csv("data/oil_data_long.csv")%>%
  group_by(family) %>%
  summarise(PERoil = mean(PERoil))%>%
  ggplot(aes(x=family, y= PERoil, fill= family))+ #community
  geom_point(size = 4, shape=21, color = "black", show.legend = F)+
  #geom_bar(stat = "identity", show.legend = F )+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "Families oil content (%)", tag = "", y= "Oil content (%)", x= "Plant family")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# summary oil content x family (some families only 1 sp)
read.csv("data/oil_data_long.csv")%>%
  group_by(community, family) %>%
  get_summary_stats(PERoil)

#### 2- ALL Oil types in percentatges ####
read.csv("data/oil_data_long.csv")%>%
  select(8:10)%>%
  group_by(oil_type)%>%
  get_summary_stats()%>%
  filter(variable == "oil_type_PER")  %>%
  #write.csv("oil_types_PER_summary.csv")%>%
  ggplot(aes(x= oil_type, y= mean, fill=oil_type)) +
  geom_point(size =5, shape = 21, color = "black", show.legend = F)+
  geom_errorbar(aes( y= mean, ymin = mean-se, ymax = mean+se),width = 0.2, size =1, color="black")+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "FA types  (%)", tag = "", y= "Oil content (%)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# barplot
read.csv("data/oil_data_long.csv")%>%
  ggplot(aes(x= oil_type, y= oil_type_PER, fill=oil_type)) +
  geom_boxplot(show.legend = F)+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "FA types  (%)", tag = "", y= "Oil content (%)", x= "FA Types")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
  
# 2.1 oil types percentatges per sp (% specific oil /total oil ID)
read.csv("data/oil_data_long.csv")%>%
  mutate(ID = paste(species, community))%>%
  group_by(Taxon, family, oil_type)%>%
  summarise(oil_type_PER = mean(oil_type_PER)) %>%
  #filter (oil_type_PER >3) %>% # if filter above 3% change geom_bar position = fill
  ggplot(aes(x=Taxon, y= oil_type_PER, fill= oil_type))+
  geom_bar(position = "stack", stat = "identity",  color = "black", show.legend = T)+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "FA types  (%)", tag = "", y= "Oil content (%)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 2.2 oil types percentatges per family (% specific oil /total oil ID)
read.csv("data/oil_data_long.csv")%>%
  group_by(family, oil_type)%>%
  summarise(oil_type_PER = mean(oil_type_PER)) %>%
 filter (oil_type_PER >2) %>% # if filter above 3% change geom_bar position = fill
  ggplot(aes(x=family, y= oil_type_PER, fill= oil_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "FA types (%)", tag = "", y= "Oil content (%)" , x= "Plant family")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

#### 3- saturated fatty acids (SFA) vs unsaturated fatty acids (UFA) ####
# per species
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/oil_type_specification.csv"), by= "oil_type")%>%
  dplyr::select(oil_type, community, Taxon, family, oil_type_PER, simple_name, saturation_type)%>%
  group_by(Taxon, family, saturation_type)%>%
  summarise(oil_type_PER = sum(oil_type_PER))%>%
  ggplot(aes(x=Taxon, y= oil_type_PER, fill= saturation_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "UFA vs SFA", tag = "E)", y= "Oil content (relative %)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))


# per family
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/oil_type_specification.csv"), by= "oil_type")%>%
  dplyr::select(oil_type, community, Taxon, family, oil_type_PER, simple_name, saturation_type)%>%
  group_by(family, saturation_type)%>%
  summarise(oil_type_PER = sum(oil_type_PER))%>%
  ggplot(aes(x=family, y= oil_type_PER, fill= saturation_type))+
  geom_bar(position = "fill", stat = "identity",  color = "black", show.legend = T)+
  coord_flip()+
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA vs SFA  (%)", tag = "E)", y= "Oil content", x= "Plant family")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# ratio UFA/SFA calculation
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/oil_type_specification.csv"), by= "oil_type")%>%
  dplyr::select(oil_type, community, Taxon, family, oil_type_PER, simple_name, saturation_type)%>%
  group_by(Taxon, family, saturation_type)%>%
  summarise(oil_type_PER = sum(oil_type_PER))%>%
  spread(saturation_type, oil_type_PER)%>%
  mutate(ratio = (UFA/SFA))%>%
  data.frame()
  group_by(family) %>%
  summarise(SFA = mean(SFA), UFA= mean(UFA), ratio = mean(ratio))

#### 4- oil content and ratio UFA/SFAvs vs predictors: seed mass,  Ecology, GDD, FDD, germ under snow ####
# seed mass calculation and included into species.csv file
read.csv("data/seed_mass.csv", sep = ";") %>%
  group_by(species, community) %>%
  summarise(seedmass= mean(mass_50))%>%
  mutate (Taxon = species) %>%
  mutate(species=make.cepnames(species))%>%
  mutate (ID= paste(Taxon, community)) -> seed_mass

# check if species names match
#setdiff(seed_mass$species,oil_per$species)
# 4.1   Oil content x seed mass
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  dplyr::select(community, species, Taxon, family, PERoil, seedmass)%>%
  group_by(species, community,  Taxon, family)%>%
  summarise(PERoil = mean(PERoil), seedmass = mean(seedmass))%>%
  #data.frame()
  ggplot(aes(x=log(seedmass), y = log(PERoil)))+ #, fill = family
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth(method= "lm")+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  theme_minimal(base_size=12) +
  labs(title= "Oil content (%) vs seed mass", tag = "", y= "Oil content (%)", x= "50 seeds mass (mg)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.2 UFA/SFA ratio x seed mass
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/oil_type_specification.csv"), by= "oil_type")%>%
  dplyr::select(oil_type, community, Taxon, family, oil_type_PER, simple_name, saturation_type)%>%
  group_by(Taxon, community, family, saturation_type)%>%
  summarise(oil_type_PER = sum(oil_type_PER))%>%
  spread(saturation_type, oil_type_PER)%>%
  mutate(ratio = (UFA/SFA))%>%
  #data.frame()
  merge(seed_mass)%>%
  dplyr::select(Taxon, family, community, ratio, seedmass)%>%
  ggplot(aes(x=log(seedmass), y = ratio))+ #, fill = family
  geom_point(shape = 21, size = 4, color="black")+
  geom_smooth(method= "lm")+
  #geom_text_repel (aes (x = seedmass, y = ratio, label = Taxon), show.legend = F, size =5, max.overlaps = 15) +
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA/SFA vs seed mass", tag = "", y= "Ratio UFA/SFA", x= "50 seeds mass (mg)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# 4.3   Oil content x ecology 
read.csv("data/data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  ggplot(aes(x=ecology, y = PERoil, fill = ecology))+
  geom_boxplot(shape = 21, size = 1, color = "black")+
  facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Oil content vs ecology", tag = "", y= "Oil content (%)", x= "Ecology")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.4   ratio UFA/SFA x ecology 
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  ggplot(aes(x=ecology, y = ratio, fill = ecology))+
  geom_boxplot(shape = 21, size = 1, color = "black")+
  #facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA/SFA vs ecology", tag = "", y= "Ratio UFA/SFA", x= "Ecology")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.5   Oil content x GDD
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family)%>%
  summarise(PERoil = mean(PERoil), GDD = mean(GDD), FDD = mean(FDD))%>%
  #data.frame()
  ggplot(aes(x=GDD, y = PERoil))+#, fill = community
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth (method = "lm")+
  #geom_text_repel(aes (x = GDD, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  #facet_grid(~community) +
  theme_minimal(base_size=12) +
  labs(title= "Oil content (%) vs GDD", tag = "", y= "Oil content (%)", x= "Growing degree days (?C)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# Oil content x FDD
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family)%>%
  summarise(PERoil = mean(PERoil), GDD = mean(GDD), FDD = mean(FDD))%>%
  #data.frame()
  ggplot(aes(x=FDD, y = PERoil))+ #, fill = community
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth (method = "lm") +
  #geom_text_repel(aes (x = FDD, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  #facet_grid(~community) +
  theme_minimal(base_size=12) +
  labs(title= "Oil content (%) vs FDD", tag = "", y= "Oil content (%)", x= "Freezing degree days (?C)")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.6  UFA/SFA ratio x GDD
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family, ecology)%>%
  summarise(ratio = mean(ratio), GDD = mean(GDD), FDD = mean(FDD))%>%
  ggplot(aes(x=GDD, y = ratio))+ #, fill = family
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth(method = "lm")+
  #geom_text_repel(aes (x = GDD, y = ratio, label = species), show.legend = F, size =5, max.overlaps = 15) +
  #facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA/SFA vs GDD", tag = "", y= "Ratio UFA/SFA", x= "GDD")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

# 4.6  UFA/SFA ratio x FDD
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family, ecology)%>%
  summarise(ratio = mean(ratio), GDD = mean(GDD), FDD = mean(FDD))%>%
  ggplot(aes(x=FDD, y = ratio))+ #, fill = community
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth(method = "lm")+
  #geom_text_repel(aes (x = FDD, y = ratio, label = species), show.legend = F, size =5, max.overlaps = 15) +
  #facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA/SFA vs FDD", tag = "", y= "Ratio UFA/SFA", x= "FDD")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.7  oil content + germ under snow
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family, ecology)%>%
  summarise(PERoil = mean(PERoil), under_snow = mean(under_snow))%>%
  ggplot(aes(x=under_snow, y = PERoil, fill = family))+
  geom_point(shape = 21, size = 4, color = "black")+
  #geom_text_repel(aes (x = under_snow, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Oil content vs under_snow", tag = "", y= "Ratio UFA/SFA", x= "Germination under snow")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))
# 4.8  ratio+ germ under snow
read.csv("data/oil_data_long.csv")%>%
  merge(read.csv("data/species.csv"))%>%
  group_by(species, community,  Taxon, family, ecology)%>%
  summarise(ratio = mean(ratio), under_snow = mean(under_snow))%>%
  ggplot(aes(x=under_snow, y = ratio))+ #, fill = family
  geom_point(shape = 21, size = 4, color = "black")+
  geom_smooth(method = "lm")+
  #geom_text_repel(aes (x = under_snow, y = ratio, label = species), show.legend = F, size =5, max.overlaps = 15) +
  #facet_grid(~community)+
  theme_minimal(base_size=12) +
  labs(title= "Ratio UFA/SFA vs under_snow", tag = "", y= "Ratio UFA/SFA", x= "Germination under snow")+
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 16),
        strip.text = element_text (size = 15),
        legend.position = "right", 
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"))

#### 5- explorative PCA with FAMe composition table #####

read.csv("data/species.csv") -> explanatory

read.csv("data/oil_per_wide.csv") %>%
  dplyr::select(code, C12.0:C24.0)%>%
  merge(explanatory, by = "code")%>%
  dplyr::select(C12.0:C24.0, PERoil, ratio)-> pca.oil


oil_all<- rda (pca.oil ~ecology+seedmass+GDD+FDD, data = explanatory, scale=TRUE)
m0 <- rda(pca.oil ~ 1, data=explanatory)

oil_all

constrained_eig <- oil_all$CCA$eig/oil_all$tot.chi*100
unconstrained_eig <- oil_all$CA$eig/oil_all$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')

RsquareAdj(oil_all)$adj.r.squared 

sel.oil_all <- ordiR2step(m0, scope = formula(oil_all), direction ="forward")

sel.oil_all$anova
vif.cca(oil_all)

ordiplot(sel.oil_all)
