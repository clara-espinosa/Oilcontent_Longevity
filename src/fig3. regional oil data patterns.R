library(tidyverse);library(rstatix);library(vegan);library(patchwork)

# FIGURE 2. REGIONAL extended oil dataset
### Panel A. Oil content differences across altitudinal pattern ###
x11()
read.csv("data/oil_regionaldata.csv", sep = ",")%>%
  convert_as_factor(family, ecology)%>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
  ggplot()+
  geom_boxplot(aes(x=ecology, y=oil_content, fill=ecology), color= "black")+
  geom_point(aes(x=ecology, y=oil_content, fill=ecology), color= "black", alpha = 0.5, position = "jitter", shape = 21, size= 4, show.legend = F)+
  scale_fill_manual (name= "Ecology distribution", values = c("orange3", "darkcyan", "orchid4"), 
                     guide = guide_legend (title.position = "top",direction = "horizontal"))+
  ggthemes::theme_tufte(base_size=12) + 
  labs( title= "A) Species oil content (%)", y= "Oil content (%)")+ #",
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 12),
        legend.position = "inside",
        legend.position.inside = c(0.5,-0.06), 
        plot.margin = unit(c(0, 0,0.2,0.1), "cm"),
        legend.title = element_text(size = 10, color = "black", hjust = 0.5),
        legend.text = element_text(size = 10, color = "black"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.box.margin = unit(c(0, 0,0,0), "cm"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"))-> fig3A;fig3A


### PANEL B) Oil content - Seed mass relationship ####
col <- c("#F2F9AA",  "#f3ec70","#f6e528","#Fad220", "#fcb622","darkgoldenrod3", "#f68b08", "#ff7125" ,"#c87107",  "#8E5005",
         "#08f9dd", "#16cacb", "#21a8be", "#2d7faf", "#275381", "#42346f",
         "#78E208", "#45A747", "#396F3E")
col_order <- c( "#8E5005","darkgoldenrod3", "gold1", "orange","darkorange2", "orangered1" ,
                "#08f9dd", "#21a8be", "#2d7faf", "#275381", "#42346f",
                "#78E208")
read.csv("data/oil_regionaldata.csv", sep = ",")%>%
  convert_as_factor(family, ecology)%>%
  mutate(ecology= fct_relevel(ecology, "Strict lowland", "Generalist", "Strict alpine"))%>%
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
  rename(seed_mass= X50seed_mass_mg)%>%
  #filter(source == "own_data")%>%
  ggplot()+
  geom_point(aes(y=seed_mass, x=oil_content, fill=order), color= "black", shape = 21, size= 4)+
  geom_smooth(aes(y=seed_mass, x=oil_content),method = "glm", color = "black", se=F)+
  scale_fill_manual (name= "Plant orders", values=col_order)+
  ggthemes::theme_tufte(base_size=12) + 
  annotate("text", label="Post mean: - 0.02\n CI: [ - 0.06 | 0.01 ]", y=300, x= 30, size = 3)+
  labs( title = "B) Oil content vs seed mass relationship", x= "Oil content (%)", y= "Seed mass (mg)")+ 
  theme(text = element_text(family = "sans"),
        plot.title = element_text (size= 12),
        legend.position = "right", 
        plot.margin = unit(c(0, 0,0.2,0.1), "cm"),
        legend.title = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 10, color = "black"),
        panel.background = element_rect(color = "black", fill = NULL),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"))-> fig3B;fig3B

# combine plots

fig3A + fig3B + plot_layout(widths= c(0.45, 0.55))-> fig3;fig3

# save plots
ggsave(filename = "fig 3. regional patterns.png", plot =fig3 , path = "results/figures", 
       device = "png",dpi = 600) #, width = 180, height = 150, units = "mm"scale = 1, 
