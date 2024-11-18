# figure ecological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr);library (viridis)

## PANEL A) OIL content and ratio vs Taxon distribution (specialist vs generalist)####
read.csv("data/species_oil.csv")%>%
  merge (oil_data)%>%
  dplyr::select(community, Taxon, family, ecology, oil.content, ratio)%>%
  group_by(community, Taxon, family, ecology)%>%
  gather(trait, value, oil.content:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(ecology = as.factor(ecology))%>%
  mutate(ecology= fct_relevel(ecology,"Strict alpine", "Generalist"))%>% 
  #data.frame()
  ggplot(aes(x=ecology, y = log(value), fill = ecology))+ #, fill = family
  geom_boxplot(size = 1, color = "black", show.legend = T)+
  scale_fill_manual (name= "Ecology", breaks = c ("Generalist", "Strict alpine"),values=c( "chocolate2", "deepskyblue3"))+
  facet_grid(~trait, scales= "free")+
  #coord_flip()+
  #geom_text_repel (aes (x = seedmass, y = oil.content, label = Taxon), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "A)", x= "Ecology",y= "log oil content (%)                                                                       Ratio UFA/SFA")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        #plot.tag.position =c(0.01,0.98), 
        #plot.margin = unit(c(0,0,0,0),'cm'),
        #legend.position = "rigth", 
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.text.x= element_blank (),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"))-> fig5A;fig5A
  

## PANEL b) OIL content and ratio vs GDD ####
read.csv("data/species_traits_summary.csv")%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  dplyr::select(community, Taxon, family, GDD, oil.content, ratio)%>%
  na.omit()%>%
  group_by(community, Taxon, family, GDD)%>%
  gather(trait, value, oil.content:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  #data.frame()
  ggplot(aes(y=value, x = GDD))+ #, fill = family
  geom_point(aes(fill = family),shape = 21, color = "black", show.legend = T, size= 4)+
   scale_fill_manual (values=col)+
  geom_smooth(method= "lm", color= "black")+
  facet_grid(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "B)")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text(size= 12),
        plot.tag.position =c(0.01,0.99), 
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "left",  
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text(12),
        axis.title.x = element_text(12),
        axis.text.x = element_text(12),
        axis.text.y = element_text(size = 10, color = "black")) -> fig5B;fig5B
## PANEL c) OIL content and ratio vs FDD ####
read.csv("data/species_traits_summary.csv")%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  dplyr::select(community, Taxon, family, FDD, oil.content, ratio)%>%
  na.omit()%>%
  group_by(community, Taxon, family, FDD)%>%
  gather(trait, value, oil.content:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  #data.frame()
  ggplot(aes(y=value, x = FDD))+ #, fill = family
  geom_point(aes(fill = family),shape = 21, color = "black", show.legend = T, size= 4)+
  scale_fill_manual (values=col)+
  geom_smooth(method= "lm", color= "black")+
  facet_grid(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "C)")+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text(size= 12),
        plot.tag.position =c(0.01,0.99), 
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "left",  
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text(12),
        axis.title.x = element_text(12),
        axis.text.x = element_text(12),
        axis.text.y = element_text(size = 10, color = "black"))-> fig5C;fig5C

# combine plots ###
ggpubr::ggarrange(fig5A, fig5B, fig5C,ncol =1, nrow= 3,common.legend = FALSE, heights = c(1,1,1),align = "v")->fig5;fig5

library(patchwork)
fig5B/fig5C + plot_layout(guides = "collect")-> fig5BC;fig5BC

fig5A/fig5BC + plot_layout(heights = c(0.5, 1))
