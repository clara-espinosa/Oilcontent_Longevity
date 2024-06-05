library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(glmmTMB);library(DHARMa)
library(phylosignal);library(phylobase);library(ape);library(tidytree)

# oil content data (own + literature)
read.csv("data/oil_fulldataset.csv")%>%
  dplyr::select(Taxon, family, ecology, seed.oil, mean.value)%>%
  filter(seed.oil=="sp.level") -> oil_fulldataset

### Phylo tree #####
# always check family names with http://www.mobot.org/MOBOT/research/APweb/
read.csv("data/oil_fulldataset.csv", sep =",") %>% 
  filter(seed.oil=="sp.level")%>%
  dplyr::select (Taxon, family) %>%
  unique %>%
  separate(Taxon, into = c("genus", "species"), sep = " ") %>%
  mutate(species = paste(genus, species),
         genus = genus,
         family = family) %>%
  arrange(species) %>%
  na.omit %>%
  select(species, genus, family)-> 
  ranks1

#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
phylo.maker(sp.list = ranks1,
            tree = GBOTB.extended, 
            nodes = nodes.info.1, 
            scenarios = "S3") ->
  tree

rm(ranks1)
x11()
plot(tree$scenario.3)

write.tree(tree$scenario.3, file = "results/tree_oil_fulldataset.tree")

# phylosignal tree graph ###
ape::read.tree("results/tree_oil_fulldataset.tree")-> tree
x11()
plot(tree)
read.csv("data/oil_fulldataset.csv", sep =",") %>% 
  filter(seed.oil=="sp.level")%>%
  rename(oil.content=mean.value)%>%
  select(Taxon, oil.content, ecology)%>%
  mutate(label= gsub(" ","_", Taxon))%>%
  remove_rownames %>% 
  column_to_rownames(var="label") %>% 
  select(oil.content)-> oil_phylo 

#https://www.francoiskeck.fr/phylosignal/demo_plots.html

obj_oil <- phylo4d(as(tree, "phylo4"), data.frame(oil_phylo), match.data=TRUE)
#mat.col <- ifelse(tdata(obj_oil , "tip") < 1,  "darkred","darkgrey")
#mat.e <- matrix(abs(rnorm(19 * 3, 0, 0.5)), ncol = 3,
# dimnames = list(tipLabels(p4d), names(tdata(p4d))))
barplot.phylo4d(obj_oil, tree.ratio = 0.2,  center=F, bar.col = rainbow(95), #bar.col = mat.col,  error.bar.sup = mat.e, 
                trait.bg.col = "white", show.box = F, grid.vertical = F,
                grid.horizontal = F, tip.labels = NULL, tree.ladderize =T) #rainbow(37)

barplot(obj_oil ,tree.type = "fan", tip.cex = 0.6, tree.open.angle = 160, trait.cex = 0.6)
dotplot.phylo4d (obj_oil , tree.type= "cladogram")
gridplot.phylo4d (obj_oil ,  tip.cex = 1.5, show.trait = T) #tree.type = "fan",
phyloSignal(obj_oil )
phyloCorrelogram (obj_oil )
  
## read species with seed mass data ####
read.csv("data/seed_mass.csv", sep = ",")%>%
  group_by(Taxon, community)%>%
  get_summary_stats(mass_50)%>%
  dplyr::select(Taxon, community, n, mean, sd)%>%
  as.data.frame()-> seedmass # 66 species
### read germ data (t50/cold_strat_germ) ####
read.csv("data/germ_t50_data.csv")%>% # so far only t50 = !!
  #na.omit()->t50  #45 species with t50 trait
  full_join(read.csv("data/germ_drivers_data.csv"))%>% # 69 species some with NAs
  na.omit()-> germ_traits 

#### Load SPECIES PREFERENCES data####
# based from 80 iButtons floristic relevés
# original scripts and data in Omaña80 and Picos Github repository (respectively)
#script name: indices and glms

# FROM ibuttons data: in picos from 02/10/2018 to 07/08/2019; in villabandin from 12/07/2021 to 29/05/2022
# FDD: sum of degrees (daily mean) when daily Tmean is below 0 ºC (in absolute values)
# GDD: sum of degrees (daily mean) when daily Tmean is above 5 ºC 
# Snowdays = sum of days with snow cover, calculated as days when Tmax <0.5 and Tmin > -0.5 ºC
# These variables (altogether with bio1, bio2 and bio7) calculated mean x day, then month and finally year (whole datset)

# From inventory data: 80 plots per each community, at each plot species coverage estimated in %
# To calculate species preferences several filters applied:
# consider only plots where each species has at least 10% coverage
# climatic variables weighted per coverage 

read.csv("data/sp_pref_picos.csv", sep =",")%>%
  rbind(read.csv("data/sp_pref_villa.csv", sep =",")) -> sp_pref # 127 species

# seed mass (log transformed) #####
oil_fulldataset%>%
  merge(seedmass)%>% #n=38
  rename(familia = family)%>%
  convert_as_factor(Taxon, community, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Lseedmass = log(mean), 
         Loil.content= log(mean.value))%>%
  na.omit()-> oil_seedmass # 
# data distribution
hist(oil_seedmass$mean) 
hist(oil_seedmass$Lseedmass)# normal distribution after log transformation
hist(oil_seedmass$mean.value)
hist(oil_seedmass$Loil.content) # almost normal distribution after log transformation


##### GLMM with dharma 
# family  Gamma(link="log")   gaussian
a <- glmmTMB(Loil.content~ Lseedmass + (1| familia) ,  family = gaussian, data= oil_seedmass) 
 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)

# Gamma(link = "log"), model plots GOOD in oil content 

oil_seedmass%>%
  dplyr::select(Taxon, familia, mean.value, mean, Lseedmass, Loil.content)%>%
  rename(oil.content=mean.value)%>%
  mutate(familia = as.factor(familia))%>%
  #mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
  #                           "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
  #                           "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
  #                            "Cyperaceae", "Juncaceae"))%>% 
  #data.frame()
  ggplot(aes(y=Lseedmass, x = Loil.content))+ #, fill = family
  geom_point(aes(fill =familia),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=rainbow(17))+
  #facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "A)", y= " log seeds mass (mg)")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,0.95), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"))

# ecological trade-offs ####
#generalist vs specialist ####
oil_fulldataset%>% # (n=95)
  rename(familia = family)%>%
  rename(oil.content = mean.value)%>%
  convert_as_factor(Taxon, familia, ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content))%>%
  na.omit()-> oil_ecology

# Gamma(link="log")
a <- glmmTMB(Loil.content  ~ ecology + (1| familia) ,  family = gaussian, data= oil_ecology) 
 
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)

oil_ecology%>%
  mutate(ecology = as.factor(ecology))%>%
  mutate(ecology= fct_relevel(ecology,"Specialist", "Generalist"))%>% 
  #data.frame()
  ggplot(aes(x=ecology, y = Loil.content, fill = ecology))+ #, fill = family
  geom_boxplot(size = 1, color = "black", show.legend = T)+
  scale_fill_manual (name= "Ecology", breaks = c ("Generalist", "Specialist"),values=c( Generalist = "chocolate2", Specialist = "deepskyblue3"))+
  #facet_grid(~trait, scales= "free")+
  coord_flip()+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "A)", x= "Ecology",y= "log oil content (%)                                                                       Ratio UFA/SFA")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,0.98), 
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "left", 
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.x= element_blank (),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"))

#species preferences ####
oil_fulldataset%>%
  merge(sp_pref)%>% #n=44
  rename(familia=family)%>%
  rename(oil.content = mean.value)%>%
  convert_as_factor(Taxon, community, familia,ecology) %>%
  mutate(ID = gsub(" ", "_", Taxon), animal = ID)%>%
  mutate(Loil.content = log(oil.content))%>%
  na.omit()-> oil_sp_pref 

# Gamma(link="log")
a <- glmmTMB(Loil.content ~ FDD+ (1| familia) ,  family = gaussian, data= oil_sp_pref)  
summary(a)
residuals <- simulateResiduals (a) ; plot(residuals)
library(viridis)
#GDD
oil_sp_pref%>%
  #data.frame()
  ggplot(aes(y=GDD, x = Loil.content, fill = GDD))+ #, fill = family
  geom_point(shape = 21, color = "black", show.legend = T, size= 4)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_viridis (option ="plasma")+
  #facet_grid(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "B)", x= "log oil content (%)                                                                       Ratio UFA/SFA")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_blank(),
        plot.tag.position =c(0.01,0.99), 
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "left",  
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black")) 
#FDD
oil_sp_pref%>%
  ggplot(aes(y=FDD, x = Loil.content, fill = FDD))+ #, fill = family
  geom_point(shape = 21, color = "black", show.legend = T, size= 4)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_viridis (option ="mako")+
  #facet_grid(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "C)",y= "", x= "Log oil content (%)                                                            Log Ratio UFA/SFA")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_blank(),
        plot.tag.position =c(0.01,0.99), 
        plot.margin = unit(c(0,0.1,0,0),'cm'),
        legend.position = "left", 
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"))
