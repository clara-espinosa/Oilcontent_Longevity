# figure biological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(ggpubr);library(ggrepel)

#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622",  "#f68b08","#ff7125" ) #, ,"#c87107",  "#8E5005"
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622", "#f68b08", "#ff7125" ,
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")

### PANEL A)  seed mass x oil content and seed mass x ratio Ufa/Sfa #####
read.csv("data/species_oil.csv")%>%
  merge(oil_data_full)%>%
  merge(seedmass)%>%
  dplyr::select(Taxon, family, mass_50, oil.content, ratio)%>%
  #mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  #data.frame()
  ggplot(aes(y=oil.content, x = mass_50))+ #, fill = family
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=col)+
  #facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "A)", title= "Oil content vs Seed mass", y= "Oil content (%)", x= "50 Seeds mass (mg)")+
  annotate("text", label="Post mean: - 0.22\n pMCMC: 0.04", x=200, y=27.5)+
  theme(text = element_text(family = "sans"),
        plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,1), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"))-> fig4A_oil;fig4A_oil

#ratio
read.csv("data/species_oil.csv")%>%
  merge(oil_data_full)%>%
  merge(seedmass)%>%
  dplyr::select(Taxon, family, mass_50, oil.content, ratio)%>%
  #mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  #data.frame()
  ggplot(aes(y=ratio, x = mass_50))+ #, fill = family
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=col)+
  #facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (x = seedmass, y = PERoil, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "", title= "Ratio UFA/SFA vs Seed mass", y= "Ratio UFA/SFA", x= "50 Seeds mass (mg)")+
  annotate("text", label ="Post mean: - 0.07\n pMCMC: 0.096", x=200,y= 11)+
  theme(text = element_text(family = "sans"),
        plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,1), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"))-> fig4A_ratio;fig4A_ratio

# combine
library(patchwork)
fig4A_oil + fig4A_ratio + 
  plot_layout(guides = 'auto')-> fig4A;fig4A


### PANEL B) viability loss x oil content and ratio Ufa/Sfa #####
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
read.csv("data/2022/genstat.csv", sep =",")%>% 
  group_by(Taxon, community)%>%
  summarise(p50=mean(p50))%>% #32
  merge(read.csv("data/species_oil.csv"), by=c("Taxon", "community")) %>% #30
  merge(oil_data_full,by=c("species", "family")) %>% #26
  dplyr::select(species, community,  family, p50, oil.content, ratio)%>%
  convert_as_factor(species, community) %>%
  ggplot(aes(x=p50, y=oil.content), color="black")+
  geom_point(aes(fill=family), size= 4, shape=21)+
  labs(title= "Oil content vs seed longevity", x= "p50 (days)", y = "Oil content (%)", tag= "B)")+
  #geom_text_repel(aes(x=PERoil, y=p50,label=species))+
  geom_smooth(aes(x=p50, y=oil.content),method="lm", color= "black")+
  scale_fill_manual (values=col)+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position =c(0.02,0.97), 
         plot.margin = unit(c(0,0,0,0),'cm'),
         axis.title.y = element_text (size=12),
         axis.text.y = element_text (size = 10, color = "black"),
         axis.title.x = element_text (size=12),
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =12, face="bold"),
         legend.key.size = unit(0.75,"line"),
         legend.margin = margin(0, 0, 0, 0),
         legend.text = element_text (size =10),
         legend.position = "none")->fig4B_oilp50; fig4B_oilp50 # legend.position = c(0.85, 0.5),
         #legend.box.background = element_rect(color = "black", size = 1)

# p50 scatterplot (genstat) ratio
read.csv("data/2022/genstat.csv", sep =",")%>% 
  group_by(Taxon, community)%>%
  summarise(p50=mean(p50))%>% #32
  merge(read.csv("data/species_oil.csv"), by=c("Taxon", "community")) %>% #30
  merge(oil_data_full,by=c("species", "family")) %>% #26
  dplyr::select(species, community,  family, p50, oil.content, ratio)%>%
  convert_as_factor(species, community) %>%
  ggplot(aes(x=p50, y=ratio), color="black")+
  geom_point(aes(fill=family), size= 4, shape=21)+
  labs(title= "Ratio vs seed longevity", x= "p50 (days)", y = "Ratio UFA/SFA", tag= "")+
  #geom_text_repel(aes(x=PERoil, y=p50,label=species))+
  geom_smooth(aes(x=p50, y=ratio),method="lm", color= "black")+
  scale_fill_manual (values=col)+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position =c(0.01,0.95), 
         plot.margin = unit(c(0,0,0,0),'cm'),
         axis.title.y = element_text (size=12),
         axis.text.y = element_text (size = 10, color = "black"),
         axis.title.x = element_text (size=12),
         axis.text.x= element_text (size = 10, color = "black"),
         strip.text = element_text( size = 18, hjust = 0),
         strip.background = element_blank(), 
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.title = element_text (size =12, face="bold"),
         legend.key.size = unit(0.75,"line"),
         legend.margin = margin(0, 0, 0, 0),
         legend.text = element_text (size =10),
         legend.position = "none")->fig4B_ratiop50; fig4B_ratiop50 

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
