# figure biological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr)

### PANEL A)  seed mass x oil content and seed mass x ratio Ufa/Sfa #####
# check normalities of data
read.csv("data/species_oil.csv")%>%
  dplyr:: select(Taxon, species, family, community, ecology, PERoil, seedmass, SFA:P50)%>%
  group_by(Taxon, species, family, community, ecology)%>%
  gather(trait, value, PERoil:P50)%>%
  ggplot()+
  geom_histogram(aes(x= log(value), fill= trait), color = "black") + # x= variable a estudiar, fill = segun que variable quieres colorear
  facet_wrap(~trait, ncol = 4, scales = "free") + # divide los gr?ficos por especie y mant?n 4 columnas
  theme_bw(base_size = 12)
# most not normal, try log transform all??

read.csv("data/species_oil.csv")%>%
  dplyr::select(community, species, Taxon, family, seedmass, PERoil, ratio)%>%
  group_by(community, species, Taxon, family, seedmass)%>%
  gather(trait, value, PERoil:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "PERoil"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  #data.frame()
  ggplot(aes(y=log(seedmass), x = log(value)))+ #, fill = family
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=rainbow(17))+
  facet_grid(~trait, scales= "free")+
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
        axis.text = element_text(size = 10, color = "black"))-> fig4A;fig4A


### PANEL b) viability loss x oil content and ratio Ufa/Sfa #####
# germination curves and p50 OIL CONTENT
# germination curves all species together
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(header) %>%
  mutate(germPER = germinated/seeds)%>%
  group_by(species, ageing )%>%
  summarise(germPER = mean(germPER), GDD = mean(GDD), PERoil = mean(PERoil))%>%
  ggplot(aes(x=ageing, y = germPER, group = species, color= PERoil))+ #GDD
  geom_smooth(method = "loess", se = FALSE, linewidth= 0.75)+
  labs (x= "Ageing time", y= "Final germination")+
  scale_color_viridis ()+
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         plot.margin = unit(c(0, 0,0,0), "cm"),
         axis.title.y = element_text (size=10),
         axis.text.y = element_text (size = 10, color = "black"),
         axis.title.x = element_text(size = 10), 
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 10, hjust = 0),
         legend.margin = margin(0, 0, 0, 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL),
         legend.title = element_text (size =10),
         legend.text = element_text (size =10),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig4B_oilcurves;fig4B_oilcurves

# p50 scatterplot (genstat) OILCONTENT
read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia,  distribution, Ki, slope, p50, GDD, PERoil)%>%
  convert_as_factor(species, code, community, familia,  distribution) %>%
  ggplot(aes(x=PERoil, y=p50, fill=PERoil), color="black")+
  geom_point(size= 4, shape=21)+
  labs(x= "Oil content (%)", y = "p50 (days)", tag= "B)")+
  #geom_text_repel(aes(x=PERoil, y=p50,label=species))+
  geom_smooth(method="lm", color= "black")+
  scale_fill_viridis (name = "Oil content (%)")+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position =c(0.01,0.95), 
         plot.margin = unit(c(0,0,0,0),'cm'),
         axis.title.y = element_text (size=12),
         axis.text.y = element_text (size = 10, color = "black"),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =12, face="bold"),
         legend.key.size = unit(0.75,"line"),
         legend.margin = margin(0, 0, 0, 0),
         legend.text = element_text (size =10),
         legend.position = "top")->fig4B_oilp50; fig4B_oilp50 # legend.position = c(0.85, 0.5),
         #legend.box.background = element_rect(color = "black", size = 1)

# combine oil content vs longevity
library(patchwork)
fig4B_oilp50+inset_element(fig4B_oilcurves, left=0.55, bottom =0.5, right=1, top=1)->fig4Boil;fig4Boil
# germination curves and p50 RATIO UFA/SFA
# germination curves all species together
read.csv("data/2022/germination22.csv", sep =";") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(header) %>%
  mutate(germPER = germinated/seeds)%>%
  group_by(species, ageing )%>%
  summarise(germPER = mean(germPER), GDD = mean(GDD), ratio = mean(ratio))%>%
  ggplot(aes(x=ageing, y = germPER, group = species, color= ratio))+ #GDD
  geom_smooth(method = "loess", se = FALSE, linewidth= 0.75)+
  labs (x= "Ageing time", y= "Final germination")+
  scale_color_viridis ()+
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         plot.margin = unit(c(0, 0,0,0), "cm"),
         axis.title.y = element_text (size=10),
         axis.text.y = element_text (size = 10, color = "black"),
         axis.title.x = element_text(size = 10), 
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 10, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL),
         legend.title = element_text (size =10),
         legend.margin = margin(0, 0, 0, 0),
         legend.text = element_text (size =10),
         legend.position = "none", # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig4B_ratiocurves;fig4B_ratiocurves

# p50 scatterplot (genstat) RATIO
read.csv("data/2022/genstat22.csv", sep =",")%>%
  left_join(header, by = c("species", "community", "code")) %>%
  dplyr::select(species, code, community, site, familia,  distribution, Ki, slope, p50, GDD, PERoil, ratio)%>%
  convert_as_factor(species, code, community, familia,  distribution) %>%
  ggplot(aes(x=ratio, y=p50, fill=PERoil), color="black")+
  geom_point(size= 4, shape=21)+
  labs(x= "Ratio UFA/SFA", y = "p50 (days)")+
  #geom_text_repel(aes(x=PERoil, y=p50,label=species))+
  geom_smooth(method="lm", color= "black")+
  scale_fill_viridis (name = "Ratio UFA/SFA")+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         #plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position = c(0,1),
         plot.margin = unit(c(0,0,0,0),'cm'),
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_blank(), 
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =12, face="bold"),
         legend.key.size = unit(0.75,"line"),
         legend.margin = margin(0, 0, 0, 0),
         legend.text = element_text (size =10),
         legend.position = "top")->fig4B_ratiop50; fig4B_ratiop50 # legend.position = c(0.85, 0.5),
#legend.box.background = element_rect(color = "black", size = 1)

# combine ratio vs longevity
library(patchwork)
fig4B_ratiop50+inset_element(fig4B_ratiocurves, left=0.55, bottom =0.5, right=1, top=1)->fig4Bratio;fig4Bratio

# combine oil + ratio
ggpubr::ggarrange(fig4Boil, fig4Bratio,ncol =2, nrow= 1,common.legend = FALSE,widths = c(1.1,1), align = "h") ->fig4B;fig4B

### PANEL C) T50 vs oil content and vs ratio #####
read.csv("data/species_oil.csv")%>%
  dplyr::select(community, species, Taxon, family, T50, PERoil, ratio)%>%
  group_by(community, species, Taxon, family, T50)%>%
  gather(trait, value, PERoil:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "PERoil"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  na.omit()%>%
  #data.frame()
  ggplot(aes(y=T50, x = log(value)))+ #, fill = family
  geom_point(aes(fill =as.factor(family)),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  facet_grid(~trait, scales= "free")+
  scale_fill_manual (values=rainbow(17))+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "C)", y= "T50 (days)")+
  theme(text = element_text(family = "sans"),
        plot.tag.position =c(0.01,0.95), 
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"))-> fig4C;fig4C

#combine all panels
ggpubr::ggarrange(fig4A, fig4B, fig4C,ncol =1, nrow= 3,common.legend = FALSE, heights = c(1,1.5,1),align = "hv")->fig4;fig4

library(patchwork)
fig4A / fig4B /fig4C +
  plot_layout(guides = 'auto')
