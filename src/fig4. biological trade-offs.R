# figure biological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(ggpubr);library(ggrepel);library(viridis)
library(patchwork)

read.csv("data/species_traits.csv")%>%
  get_summary_stats(EHS_mean)

col_order <- c( "#8E5005","darkgoldenrod3", "gold1", "orange","darkorange2", "orangered1" ,
                "skyblue1","deepskyblue",  "#2d7faf", "#275381", "#42346f",
                "#78E208")

### PANEL A)  seed mass x oil content and seed mass x ratio Ufa/Sfa #####
# with log seed mass looks better
read.csv("data/species_traits.csv")%>% # from script 2 header data handling
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, family,order,oil.content, ratio, mass_50 )%>%
  na.omit()%>%
  group_by(Taxon, family, order)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
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
  ggplot(aes(x=oil.content, y = log(mass_50)))+ #, fill = family   
  geom_point(aes(fill =order),shape = 21, size = 4, color = "black")+
  geom_smooth(method= "lm", color= "black", se=F)+
  scale_fill_manual (values=col_order)+
  scale_y_continuous (limits = c(-0.75, 6))+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag= "(a)", title= "Seed mass (n=47)", x= "Oil content (%)", y= "Seed mass (log)")+ #tag = "A)",
  annotate("text", label="Post. mean: - 0.02\n CI: [ - 0.06 | 0.02 ]", x=29, y=5.5, size=3)+ # with log x = 5 , without log x= 200
  theme(text = element_text(family = "sans"),
        plot.title= element_text( size= 12, face = "bold"), #hjust = 0.5,
        plot.tag = element_text(face="bold"),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        legend.position = "none", 
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"))-> fig4A_oil;fig4A_oil

#ratio
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, family, order, oil.content, ratio, mass_50 )%>%
  na.omit()%>%
  group_by(Taxon, family, order)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales",
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", 
                             "Saxifragales", "Poales"))%>% 
  ggplot(aes(x=ratio, y = log(mass_50)))+ #, fill = family
  geom_point(aes(fill =order),shape = 21, size = 4, color = "black")+
  geom_smooth(method= "lm", color= "black", se=F)+
  scale_fill_manual (values=col_order)+
  scale_y_continuous (limits = c(-0.75, 6))+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "", title= "", x= "Ratio UFA/SFA")+ #, y= "50 Seeds mass (log)"
  annotate("text", label ="Post mean: - 0.01\n CI: [ - 0.11 | 0.09 ]", x=17.5,y= 5.5, size=3)+
  theme(text = element_text(family = "sans"),
        plot.margin = unit(c(0, 0,0,0), "cm"),
        legend.position = "none", 
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black"))-> fig4A_ratio;fig4A_ratio

# combine
fig4A_oil + fig4A_ratio + 
  plot_layout(guides = 'auto')-> fig4A;fig4A


### PANEL B) viability loss x oil content and ratio UFA/SFA #####
# p50 scatterplot (genstat) OILCONTENT
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, community , order, family, ecology, oil.content, ratio, p50)%>% 
  group_by(Taxon, community, family, order)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), p50 = mean(p50))%>% #N=29
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales", #
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", #
                             "Saxifragales", "Poales"))%>%
  ggplot(aes(y=p50, x=oil.content), color="black")+
  geom_point(aes(fill=order), size= 4, shape=21)+
  labs(tag= "(b)",title= "Seed longevity (n=33)", y= "p50 (days)", x = "Oil content (%)")+ #, tag= "B)"
  geom_smooth(aes(y=p50, x=oil.content),method="lm", color= "black", se = F)+
  annotate("text", label ="Post mean: - 0.85\n CI: [ - 1.31 | - 0.39 ]", y=45,x= 29, size=3)+
  scale_fill_manual (values=col_order)+ #direction =-1
  scale_y_continuous (limits = c(0,50))+
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 12), #hjust = 0.5,
         plot.tag = element_text(face="bold"),
         plot.margin = unit(c(0,0,0.5,0),'cm'),
         axis.title = element_text (size=10),
         axis.text = element_text (size = 9, color = "black"),
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.position = "none")->fig4B_oil; fig4B_oil # legend.position = c(0.85, 0.5),
         #legend.box.background = element_rect(color = "black", size = 1)


# p50 scatterplot (genstat) ratio
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(Taxon, community , family,order,  ecology, oil.content, ratio, p50)%>% 
  group_by(Taxon, community, family, order)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), p50 = mean(p50))%>% 
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales", #
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", #
                             "Saxifragales", "Poales"))%>%
  convert_as_factor(Taxon, community) %>%
  ggplot(aes(y=p50, x=ratio), color="black")+
  geom_point(aes(fill=order), size= 4, shape=21)+
  labs(title= "", y= "p50 (days)", x = "Ratio UFA/SFA", tag= "")+
  geom_smooth(aes(y=p50, x=ratio),method="lm", color= "black", se = F)+
  annotate("text", label ="Post mean: - 0.73\n CI: [ -2.29 | 1.15 ]", y=45,x= 17.5, size=3)+
  scale_y_continuous (limits = c(0,50))+
  scale_fill_manual (values=col_order)+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.margin = unit(c(0,0,0.5,0),'cm'),
         axis.title.y = element_blank (),
         axis.text.y = element_blank (),
         axis.title.x = element_text (size=10),
         axis.text.x= element_text (size = 9, color = "black"),
         panel.background = element_rect(color = "black", fill = NULL), 
         legend.position = "none")->fig4B_ratio; fig4B_ratio 

# combine 
fig4B_oil + fig4B_ratio+ 
  plot_layout(guides = 'auto')->fig4B;fig4B

### PANEL C) EHS vs oil content and vs ratio #####
# oil content
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(community, Taxon,order, family, EHS_mean, oil.content, ratio)%>%
  group_by(Taxon, community, family, order)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), EHS = mean(EHS_mean))%>% 
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales", #
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", #
                             "Saxifragales", "Poales"))%>%
  ggplot(aes(y=EHS, x = oil.content))+ #, fill = family
  geom_point(aes(fill =as.factor(order)),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black", se= F)+
  scale_fill_manual (values=col_order,guide = guide_legend(nrow = 3) )+
  annotate("text", label ="Post mean: 0.02\n CI: [ - 0.01 | 0.05 ]", x=29,y= 1150, size=3)+
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "(c)", title= "Germination timing (n=34)", x= "Oil content (%)", y= "EHS (ºC)")+
  theme(text = element_text(family = "sans"),
        plot.title= element_text( size= 12, face = "bold"),
        plot.tag = element_text(face="bold"),
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 9, color = "black"),
        legend.key.size = unit(1,"line"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9, color = "black"))-> fig4C_oil;fig4C_oil
# ratio
read.csv("data/species_traits.csv")%>%
  merge(read.csv("data/species_header.csv"), by = c ("Taxon", "community"))%>%
  dplyr::select(community, Taxon,order, family, EHS_mean, oil.content, ratio)%>%
  group_by(Taxon, community, family, order)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), EHS = mean(EHS_mean))%>% 
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Asteraceae","Campanulaceae",  "Apiaceae","Lamiaceae", "Orobanchaceae" ,
                              "Plantaginaceae","Gentianaceae" ,"Primulaceae","Caryophyllaceae", "Plumbaginaceae",  
                              "Fabaceae","Salicaceae","Brassicaceae", "Cistaceae", 
                              "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  mutate(order = as.factor(order))%>%
  mutate(order = fct_relevel(order,"Asterales", "Apiales","Lamiales", #
                             "Gentianales" ,"Ericales","Caryophyllales",  
                             "Fabales","Malpighiales","Brassicales", "Malvales", #
                             "Saxifragales", "Poales"))%>%
  ggplot(aes(y=EHS, x = ratio))+ #, fill = family
  geom_point(aes(fill =as.factor(order)),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black", se= F)+
  #facet_grid(~trait, scales= "free")+
  scale_fill_manual (values=col_order)+
  annotate("text", label ="Post mean: - 0.04\n CI: [ - 0.03 | 0.14 ]", x=17.5,y= 1150, size=3)+
  labs(title = "")+
  ggthemes::theme_tufte(base_size=12) + 
  labs(x= "Ratio UFA/SFA" )+
  theme(text = element_text(family = "sans"),
        plot.title= element_text( size= 12, face = "bold"),
        plot.margin = unit(c(0,0,0,0),'cm'),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 9, color = "black"),
        legend.key.size = unit(1,"line"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(size = 9, color = "black"))-> fig4C_ratio;fig4C_ratio

#combine
(fig4C_oil + fig4C_ratio )+ 
  plot_layout(guides = "collect")& theme(legend.position = "bottom") -> fig4C;fig4C # 

x11()
fig4A / fig4B /fig4C +
  plot_layout(heights = c(1,1,1))->fig4;fig4

# save plots
ggsave(filename = "fig 4 biological tradeoffs.png", plot =fig4 , path = "results/figures", 
       device = "png", dpi = 600, width = 180, units = "mm")
