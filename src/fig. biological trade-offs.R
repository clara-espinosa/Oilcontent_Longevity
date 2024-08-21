# figure biological tradeoffs
library(tidyverse);library(readxl);library(rstatix)
library(vegan);library(ggpubr);library(ggrepel);library(viridis)

#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622",  "#f68b08","#ff7125" ) #, ,"#c87107",  "#8E5005"
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622", "#f68b08", "#ff7125" ,
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")

### PANEL A)  seed mass x oil content and seed mass x ratio Ufa/Sfa #####
# with log seed mass looks better
read.csv("data/species_traits_summary.csv")%>%
  dplyr::select(Taxon, family, oil.content, ratio, mass_50 )%>%
  na.omit()%>%
  group_by(Taxon, family)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
  #print(n=34)
  #mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  #print(n=34)
  #data.frame()
  ggplot(aes(x=oil.content, y = log(mass_50)))+ #, fill = family   
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=col)+
  scale_y_continuous (limits = c(-0.75, 6))+
  #facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (y = log(mass_50), x = oil.content, label = Taxon), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs( title= "A) Seed mass", x= "Oil content (%)", y= "50 Seeds mass (log)")+ #tag = "A)",
  annotate("text", label="Post mean: - 0.22\n pMCMC: 0.04", x=5.5, y=-0.5)+ # with log x = 5 , without log x= 200
  theme(text = element_text(family = "sans"),
        plot.title= element_text( size= 14, face = "bold"), #hjust = 0.5,
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,1), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"))-> fig4A_oil;fig4A_oil

#ratio
read.csv("data/species_traits_summary.csv")%>%
  dplyr::select(Taxon, family, oil.content, ratio, mass_50 )%>%
  na.omit()%>%
  group_by(Taxon, family)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
  #mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  #data.frame()
  ggplot(aes(x=ratio, y = log(mass_50)))+ #, fill = family
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=col)+
  scale_y_continuous (limits = c(-0.75, 6))+
  #facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (x = seedmass, y = oil.content, label = species), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs(tag = "", title= "", x= "Ratio UFA/SFA")+ #, y= "50 Seeds mass (log)"
  annotate("text", label ="Post mean: - 0.07\n pMCMC: 0.096", x=4,y= -0.5)+
  theme(text = element_text(family = "sans"),
        plot.title= element_text(hjust = 0.5, size= 14, face = "bold"),
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.01,1), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"))-> fig4A_ratio;fig4A_ratio

# combine
library(patchwork)
fig4A_oil + fig4A_ratio + 
  plot_layout(guides = 'auto')-> fig4A;fig4A

# with facet wrap
read.csv("data/species_traits_summary.csv")%>%
  dplyr::select(Taxon, family, oil.content, ratio, mass_50 )%>%
  na.omit()%>%
  group_by(Taxon, family)%>%
  summarise(mass_50 = mean(mass_50), oil.content= mean(oil.content), ratio=mean(ratio))%>%
  gather(trait, oil,oil.content:ratio)%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", 
                              "Fabaceae", "Salicaceae","Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>%
  #print(n=34)
  #data.frame()
  ggplot(aes(x=oil, y = log(mass_50)))+ #, fill = family   
  geom_point(aes(fill =family),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  scale_fill_manual (values=col)+
  scale_y_continuous (limits = c(-0.75, 6))+
  facet_grid(~trait, scales= "free")+
  #geom_text_repel (aes (y = log(mass_50), x = oil.content, label = Taxon), show.legend = F, size =5, max.overlaps = 15) +
  ggthemes::theme_tufte(base_size=12) + 
  labs( tag= "A) Seed mass", x= "Oil content (%)", y= "50 Seeds mass (log)")+ #tag = "A)",
  #annotate("text", label="Post mean: - 0.22\n pMCMC: 0.04", x=5.5, y=-0.5)+ # with log x = 5 , without log x= 200
  theme(text = element_text(family = "sans"),
        plot.title= element_text( size= 14, face = "bold"), #hjust = 0.5,
        strip.text = element_text (family = "sans",size= 12, face = "bold"),
        plot.tag.position =c(0.08,1), 
        legend.position = "none", 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"))  #-> fig4A_grid;fig4A_grid

### PANEL B) viability loss x oil content and ratio Ufa/Sfa #####
# germination curves and p50 OIL CONTENT
# germination curves all species together
read.csv("data/longevity/germination.csv", sep =",") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(read.csv("data/longevity/species.csv", sep =","), by= "code") %>%
  group_by(community,family,Taxon, ageing,year)%>%
  summarise(seeds = sum(seeds), germinated = sum(germinated))%>%
  merge(oil_data, by=c("Taxon", "community", "family"))%>% 
  mutate(germPER = germinated/seeds)%>%
  group_by(Taxon, ageing )%>%
  summarise(germPER = mean(germPER),  oil.content = mean(oil.content), ratio = mean(ratio))%>%
  ggplot(aes(x=ageing, y = germPER, group = Taxon, color= oil.content))+ #GDD
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
         legend.position = "none", 
         # legend.position = c(0.85, 0.5),
         legend.box.background = element_rect(color = "black", size = 2))-> fig4B_oilcurves;fig4B_oilcurves

# p50 scatterplot (genstat) OILCONTENT
poales <- c("#78E208", "#45A747", "#396F3E")
rosids2 <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#fcb622",  "#f68b08","#ff7125" ) #,"#f4bc50"= "Apiaceae",, ,"#c87107",  "#8E5005"
asterids2 <- c("#08f9dd", "#16cacb" ,"#21a8be",   "#275381","#42346f" ) # "#16cacb" ="Brassicaceae","#21a8be" ="Fabaceae","#2d7faf" = "Salicaceae","#42346f"="Saxifragaceae",,,,
col2 <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#fcb622", "#f68b08", "#ff7125" ,
          "#08f9dd", "#16cacb" ,"#21a8be",   "#275381","#42346f", 
         "#78E208", "#45A747", "#396F3E")

read.csv("data/longevity/genstat.csv")%>%
  merge(oil_data, by=c("Taxon", "community"))%>% 
  dplyr::select(Taxon, community , family, ecology, oil.content, ratio, p50)%>% 
  na.omit()%>% 
  group_by(Taxon, community, family)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), p50 = mean(p50))%>% #N=29
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                               "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", "Fabaceae",
                               "Crassulaceae",  "Saxifragaceae","Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  convert_as_factor(Taxon, community) %>%
  ggplot(aes(y=p50, x=oil.content), color="black")+
  geom_point(aes(fill=family), size= 4, shape=21)+
  labs(title= "B) Seed longevity", y= "p50 (days)", x = "Oil content (%)")+ #, tag= "B)"
  #geom_text_repel(aes(x=oil.content, y=p50,label=Taxon))+
  geom_smooth(aes(y=p50, x=oil.content),method="lm", color= "black", se = F)+
  annotate("text", label ="Post mean: - 0.38\n pMCMC: 0.05", y=2.5,x= 5)+
  scale_fill_manual (values=col2)+ #direction =-1
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


# combine oil content vs longevity
# combine
library(patchwork)
fig4B_oilp50+inset_element(fig4B_oilcurves, left=0.55, bottom =0.5, right=1, top=1)->fig4Boil;fig4Boil
#fig4B_oilp50 + fig4B_ratiop50  + 
  #plot_layout(guides = 'auto')-> fig4B_p50;fig4B_p50


# germination curves and p50 RATIO UFA/SFA
# germination curves all species together
read.csv("data/longevity/germination.csv", sep =",") %>%
  gather(scores, germinated, D7:D28) %>%
  group_by(code, ageing, seeds) %>%
  summarise(germinated = sum(germinated, na.rm = TRUE))%>%
  merge(read.csv("data/longevity/species.csv", sep =","), by= "code") %>%
  group_by(community,family,Taxon, ageing,year)%>%
  summarise(seeds = sum(seeds), germinated = sum(germinated))%>%
  merge(oil_data, by=c("Taxon", "community", "family"))%>% 
  mutate(germPER = germinated/seeds)%>%
  group_by(Taxon, ageing )%>%
  summarise(germPER = mean(germPER),  oil.content = mean(oil.content), ratio = mean(ratio))%>%
  ggplot(aes(x=ageing, y = germPER, group = Taxon, color= ratio))+ #GDD
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

# p50 scatterplot (genstat) ratio
read.csv("data/longevity/genstat.csv")%>%
  merge(oil_data, by=c("Taxon", "community"))%>% 
  dplyr::select(Taxon, community , family, ecology, oil.content, ratio, p50)%>% 
  na.omit()%>% 
  group_by(Taxon, community, family)%>%
  summarise(oil.content = mean(oil.content), ratio = mean(ratio), p50 = mean(p50))%>% #N=29
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Plantaginaceae","Lamiaceae", "Primulaceae", "Cistaceae", "Brassicaceae", "Fabaceae",
                              "Crassulaceae",  "Saxifragaceae","Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% 
  convert_as_factor(Taxon, community) %>%
  ggplot(aes(y=p50, x=ratio), color="black")+
  geom_point(aes(fill=family), size= 4, shape=21)+
  labs(title= "", y= "p50 (days)", x = "Ratio UFA/SFA", tag= "")+
  #geom_text_repel(aes(x=ratio, y=p50,label=Taxon))+
  geom_smooth(aes(y=p50, x=ratio),method="lm", color= "black", se = F)+
  annotate("text", label ="Post mean: - 0.12\n pMCMC: 0.13", y=5,x= 5)+
  scale_fill_manual (values=col2)+ #direction =-1
  ggthemes::theme_tufte(base_size=12) + 
  theme (text = element_text(family = "sans"),
         plot.title = element_text (face = "bold",size = 20), #hjust = 0.5,
         plot.tag.position =c(0.01,0.95), 
         plot.margin = unit(c(0,0,0,0),'cm'),
         axis.title.y = element_blank (),
         axis.text.y = element_blank (),
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

# combine ratio vs longevity
library(patchwork)
fig4B_ratiop50+inset_element(fig4B_ratiocurves, left=0.55, bottom =0.5, right=1, top=1)->fig4Bratio;fig4Bratio

# combine oil + ratio
ggpubr::ggarrange(fig4Boil, fig4Bratio,ncol =2, nrow= 1,common.legend = FALSE,widths = c(1.1,1), align = "h") ->fig4B;fig4B

### PANEL C) T50 vs oil content and vs ratio #####
#order like family fct relevel
poales <- c("#78E208", "#45A747", "#396F3E")
rosids <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622",  "#f68b08","#ff7125" ) #, ,"#c87107",  "#8E5005"
asterids <- c("#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f")
col3 <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#f4bc50","#fcb622", "#f68b08",
         "#08f9dd", "#16cacb", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")
read.csv("data/species_traits_summary.csv")%>%
  dplyr::select(community, Taxon, family, mean_T50, oil.content, ratio)%>%
  group_by(community, Taxon, family, mean_T50)%>%
  gather(trait, value, oil.content:ratio)%>%
  mutate(trait = as.factor(trait))%>%
  mutate(trait = recode(trait, "oil.content"= "Oil content (%)", "ratio"= "Ratio UFA/SFA"))%>%
  mutate(family = as.factor(family))%>%
  mutate(family = fct_relevel(family,"Caryophyllaceae", "Plumbaginaceae", "Asteraceae", "Campanulaceae",
                              "Apiaceae", "Plantaginaceae","Lamiaceae",  "Cistaceae", "Brassicaceae", 
                               "Crassulaceae", "Saxifragaceae", "Poaceae", 
                              "Cyperaceae", "Juncaceae"))%>% #"Primulaceae", "Fabaceae","Salicaceae",
  na.omit()%>%
  #data.frame()
  ggplot(aes(y=mean_T50, x = log(value)))+ #, fill = family
  geom_point(aes(fill =as.factor(family)),shape = 21, size = 4, color = "black", show.legend = T)+
  geom_smooth(method= "lm", color= "black")+
  facet_grid(~trait, scales= "free")+
  scale_fill_manual (values=col3)+
  geom_text_repel (aes (x = log(value), y = mean_T50, label = Taxon), show.legend = F, max.overlaps = 15) +
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
x11()
library(patchwork)
fig4A / fig4B /fig4C +
  plot_layout(guides = 'auto')
