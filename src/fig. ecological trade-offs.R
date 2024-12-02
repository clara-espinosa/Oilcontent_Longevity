# figure ecological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(ggrepel):library(vegan);library(ggpubr);library (viridis)

read.csv("data/species_traits_summary.csv")-> summary_table
summary_table%>%
  get_summary_stats(Snw)

## PANEL A) OIL content and ratio vs GDD ####
sig.label <- as.data.frame(cbind(trait = c("Oil content (%)", "Ratio UFA/SFA"),
                                 x= c(2100, 2100), y= c(30, 18), 
                                 label = c("pMCMC: 0.32", "pMCMC: 0.54")))
sig.label%>%
  mutate(x= as.numeric(x),
         y= as.numeric(y))->sig.label
str(sig.label)

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
  geom_point(aes(fill = family),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col)+
  geom_smooth(method= "lm", color= "black", se=F)+
  facet_wrap(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "A)", x= "Growing Degree Days (GDD)")+
  geom_text (data= sig.label, aes(x= x, y= y, label = label), size= 4.5)+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text(size= 12, face = "bold"),
        plot.tag.position =c(0.01,0.99), 
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.title.x = element_text(12),
        axis.text.x = element_text(12),
        axis.text.y = element_text(size = 10, color = "black")) -> fig5A;fig5A 

## PANEL b) OIL content and ratio vs FDD ####

sig.label2 <- as.data.frame(cbind(trait = c("Oil content (%)", "Ratio UFA/SFA"),
                                 x= c(150, 150), y= c(30, 18), 
                                 label = c("pMCMC: 0.5", "pMCMC: 0.68")))
sig.label2%>%
  mutate(x= as.numeric(x),
         y= as.numeric(y))->sig.label2
str(sig.label2)

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
  geom_point(aes(fill = family),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col)+
  geom_smooth(method= "lm", color= "black", se=F)+
  facet_wrap(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "B)", x = "Freezing Degree Days (FDD)")+
  geom_text (data= sig.label2, aes(x= x, y= y, label = label), size= 4.5)+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text(size= 12, face= "bold" ),
        plot.tag.position =c(0.01,0.99), 
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        legend.title = element_text(size = 12, color = "black", face="bold"),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.title.x = element_text(12),
        axis.text.x = element_text(12),
        axis.text.y = element_text(size = 10, color = "black"))-> fig5B;fig5B
## PANEL c) OIL content and ratio vs Snow ####
sig.label3 <- as.data.frame(cbind(trait = c("Oil content (%)", "Ratio UFA/SFA"),
                                  x= c(140, 140), y= c(30, 18), 
                                  label = c("pMCMC: 0.17", "pMCMC: 0.09")))
sig.label3%>%
  mutate(x= as.numeric(x),
         y= as.numeric(y))->sig.label3
str(sig.label3)
read.csv("data/species_traits_summary.csv")%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  dplyr::select(community, Taxon, family, Snw, oil.content, ratio)%>%
  na.omit()%>%
  group_by(community, Taxon, family, Snw)%>%
  gather(trait, value, oil.content:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  #data.frame()
  ggplot(aes(y=value, x = Snw))+ #, fill = family
  geom_point(aes(fill = family),shape = 21, color = "black", show.legend = T, size= 4)+
  scale_fill_manual (values=col, guide=guide_legend(nrow = 3))+
  geom_smooth(method= "lm", color= "black", se=F)+
  facet_wrap(~trait, scales= "free")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "C)", x= "Snow days")+
  geom_text (data= sig.label3, aes(x= x, y= y, label = label), size= 4.5)+
  theme(text = element_text(family = "sans"),
        #plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text(size= 12, face= "bold"),
        plot.tag.position =c(0.01,0.99), 
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "bottom",  
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        legend.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_blank(),
        axis.title.x = element_text(12),
        axis.text.x = element_text(12),
        axis.text.y = element_text(size = 10, color = "black"))-> fig5C;fig5C

# combine plots ##

library(patchwork)
fig5A/fig5B/fig5C + plot_layout(heights = c(1, 1, 1), guides = 'collect')&theme(legend.position = 'bottom')-> fig5;fig5

# save plot
ggsave(filename = "fig 5 ecological tradeoffs.png", plot =fig5 , path = "results/figures", 
       device = "png", dpi = 600)
