# figure ecological tradeoffs
library(tidyverse);library(readxl);library(rstatix); library(patchwork)
library(ggrepel):library(vegan);library(ggpubr);library (viridis)

col_order <- c( "#8E5005","darkgoldenrod3", "gold1", "orange","darkorange2", "orangered1" ,
                "skyblue1","deepskyblue",  "#2d7faf", "#275381", "#42346f",
                "#78E208")
## PANEL A) OIL content and ratio vs GDD ####
# oil content
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=oil.content, x = GDD))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "(a)", x= "Growing Degree Days (GDD, ºC)", y= "Oil content (%)")+
  annotate("text", label="Post. mean: - 0.0002\n CI: [ - 0.0006 | 0.0002 ]", x=1950, y=32, size=3)+
  theme(text = element_text(family = "sans"),
        plot.tag.position =c(0.01,0.99), 
        plot.tag = element_text(face="bold"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black")) -> fig5A_oil;fig5A_oil

# ratio
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=ratio, x = GDD))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs(x= "Growing Degree Days (GDD, ºC)", y= "UFA/SFA ratio")+
  annotate("text", label="Post. mean: - 0.00006\n CI: [ - 0.0002 | 0.0002 ]", x=1950, y=18.8, size=3)+
  theme(text = element_text(family = "sans"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black")) -> fig5A_ratio;fig5A_ratio

# combine
fig5A_oil + fig5A_ratio + 
  plot_layout(guides = 'auto', axis_titles = "collect_x")-> fig5A;fig5A

## PANEL b) OIL content and ratio vs FDD ####
# oil content
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=oil.content, x = FDD))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "(b)", x= "Freezing Degree Days (FDD, ºC)", y= "Oil content (%)")+
  annotate("text", label="Post. mean: 0.002\n CI: [ - 0.003 | 0.008 ]", x=138, y=29, size=3)+
  theme(text = element_text(family = "sans"),
        plot.tag.position =c(0.01,0.99), 
        plot.tag = element_text(face="bold"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black"))-> fig5B_oil;fig5B_oil

# ratio 
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=ratio, x = FDD))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = F, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs( x= "Freezing Degree Days (FDD, ºC)", y= "UFA/SFA ratio")+
  annotate("text", label="Post. mean: 0.0008\n CI: [ - 0.002 | 0.004 ]", x=138, y=17, size=3)+
  theme(text = element_text(family = "sans"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none",  
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black"))-> fig5B_ratio;fig5B_ratio

# combine
fig5B_oil + fig5B_ratio + 
  plot_layout(guides = 'auto', axis_titles = "collect_x")-> fig5B;fig5B
## PANEL c) OIL content and ratio vs Snow ####
# oil content
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=oil.content, x = Snw))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = T, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "(c)", x= "Snow cover (days)", y= "Oil content (%)")+
  annotate("text", label="Post. mean: 0.003\n CI: [ - 0.001 | 0.007 ]", x=125, y=28, size=3)+
  theme(text = element_text(family = "sans"),
        plot.tag.position =c(0.01,0.99), 
        plot.tag = element_text(face="bold"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none", 
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black"))-> fig5C_oil;fig5C_oil

# ratio
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  dplyr::select(community, Taxon, order, family, GDD, FDD, Snw, oil.content, ratio)%>%
  na.omit()%>%
  ggplot(aes(y=ratio, x = Snw))+ #, fill = family
  geom_point(aes(fill = order),shape = 21, color = "black", show.legend = T, size= 4)+
  scale_fill_manual (values=col_order)+
  geom_smooth(method= "lm", color= "black", se=F)+
  ggthemes::theme_tufte(base_size=12) + 
  labs( x= "Snow cover (days)", y= "UFA/SFA ratio")+
  annotate("text", label="Post. mean: 0.002\n CI: [ - 0.0003 | 0.004 ]", x=125, y=16.5, size=3)+
  theme(text = element_text(family = "sans"),
        #plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "none", 
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(9),
        axis.text = element_text(size = 9, color = "black"))-> fig5C_ratio;fig5C_ratio

# combine

(fig5C_oil + fig5C_ratio) + 
  
  plot_layout(guides = 'collect', axis_titles = "collect_x")& 
  theme(legend.position = 'bottom', legend.title = element_blank())-> fig5C;fig5C

# FINAL combine plots ##
fig5A/fig5B/fig5C + 
  plot_annotation(title= "Ecological drivers")+
  plot_layout(heights = c(1, 1, 1), guides = 'collect')&
  theme(plot.title = element_text(size=14, hjust = 0.5))-> fig5;fig5

# save plot
ggsave(filename = "fig 5 ecological tradeoffs.png", plot =fig5 , path = "results/figures", 
       device = "png", dpi = 600, width = 180, units = "mm")
