# figure ecological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr);library (viridis)

## PANEL A) OIL content and ratio vs species distribution (specialist vs generalist)####
read.csv("data/species_oil.csv")%>%
  dplyr::select(community, species, ecology, PERoil, ratio)%>%
  group_by(community, species, ecology)%>%
  gather(trait, value, PERoil:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "PERoil"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(ecology = as.factor(ecology))%>%
  mutate(ecology= fct_relevel(ecology,"Specialist", "Generalist"))%>% 
  #data.frame()
  ggplot(aes(x=ecology, y = log(value), fill = ecology))+ #, fill = family
  geom_boxplot(size = 1, color = "black", show.legend = T)+
  scale_fill_manual (name= "Ecology", breaks = c ("Generalist", "Specialist"),values=c( Generalist = "chocolate2", Specialist = "deepskyblue3"))+
  facet_grid(~trait, scales= "free")+
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
        axis.text.y = element_text(size = 10, color = "black"))-> fig5A;fig5A
  

## PANEL b) OIL content and ratio vs GDD ####
read.csv("data/species_oil.csv")%>%
  dplyr::select(community, species, family, GDD, PERoil, ratio)%>%
  group_by(community, species, family, GDD)%>%
  gather(trait, value, PERoil:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "PERoil"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  #data.frame()
  ggplot(aes(y=GDD, x = log(value), fill = GDD))+ #, fill = family
  geom_point(shape = 21, color = "black", show.legend = T, size= 4)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_viridis (option ="plasma")+
  facet_grid(~trait, scales= "free")+
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
        axis.text.y = element_text(size = 10, color = "black")) -> fig5B;fig5B
## PANEL c) OIL content and ratio vs FDD ####
read.csv("data/species_oil.csv")%>%
  dplyr::select(community, species, family, FDD, PERoil, ratio)%>%
  group_by(community, species, family, FDD)%>%
  gather(trait, value, PERoil:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "PERoil"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  #data.frame()
  ggplot(aes(y=FDD, x = log(value), fill = FDD))+ #, fill = family
  geom_point(shape = 21, color = "black", show.legend = T, size= 4)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_viridis (option ="mako")+
  facet_grid(~trait, scales= "free")+
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
        axis.text = element_text(size = 10, color = "black"))-> fig5C;fig5C

# combine plots ###
ggpubr::ggarrange(fig5A, fig5B, fig5C,ncol =1, nrow= 3,common.legend = FALSE, heights = c(1,1,1),align = "v")->fig5;fig5
