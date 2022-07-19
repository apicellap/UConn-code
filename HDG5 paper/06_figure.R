
#download data from gene expression survey of hdg5 variants: 
urlfile4a="https://raw.githubusercontent.com/apicellap/data/main/hdg5_variants_fig6.csv"
hdg5vars<-read.csv(url(urlfile4a)) 
head(hdg5vars)


View(hdg5vars)
ctl <- expression(paste(bold("Control")))
gene1 <- expression(paste(bolditalic("HDG5-7")))
gene2 <- expression(paste(bolditalic("HDG5-8")))
gene3 <- expression(paste(bolditalic("HDG5-9")))

library(ggplot2)

hdg5.plot <-ggplot(hdg5vars, aes(x = gene, y = exp)) + 
  geom_bar(stat = 'identity', #instructs R not to make any additional calculations on the data 
           fill = "#82807E", #gray color hexcode for the bar fill
           color = "black", #color of bar border
           size =1) + #thickness of bar border
  geom_errorbar(aes(x = gene, ymin = exp - sd_m, ymax = exp + sd_p, width=0.1), size = 1) + #error bar proportions
  xlab("") + ylab("Relative expression")  +
  scale_x_discrete(limits = c("control", "hdg5-7", "hdg5-8", "hdg5-9"), #provide names of the hdg5 variants in the dataframe
                   labels = c(ctl, gene1, gene2, gene3)) + #provide object names (created above) that correspond to the dataframe info
  scale_y_continuous(
    limits = c(0,24), #set y axis range 
    breaks= seq(0,24, by = 4)) + 
  theme(
    axis.line.y.left = element_line(color = "black"), #color of left axis line
    axis.title.y.left = element_text(color="black", face="bold", size=18), #color, face, and size of y axis title
    axis.text.y.left = element_text(color="black", face="bold", size =18), #color, face, and size of y axis values
    axis.text.x = element_text(color="black", face="italic", size = 18),  #color, face, and size of y axis values
    axis.ticks = element_line(size = 1), #size of tick on x and y axes
    axis.ticks.length = unit(-0.25, "cm"), #makes ticks go inside the plot 
    axis.title.x=element_text(color="black",size=18), #color, size of x axis title (not included)
    plot.title = element_text(face="bold", size=19, hjust = 0.5), 
    plot.caption = element_text(size=12.5, colour ="red"),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x = 3, y = 22, label = "***", size = 8, color = "black") + 
  annotate("text", x = 4, y = 6, label = "**", size = 8, color = "black") 


hdg5.plot 

urlfile1b="https://raw.githubusercontent.com/apicellap/data/main/tri_count.csv"
tri_ct<-read.csv(url(urlfile1b)) 

library(dplyr)
library(ggpubr)

#create summary statistics table (a new dataframe called 'tri_mean')
tri_mean<-summarise(
  group_by(tri_ct, tissue, tri_type), 
  obs = n(), #number of observations (biological reps) per group
  Mcount = mean(count), 
  SEcount =(sd(count))/sqrt(obs)) #multiple statistics can be calculated within summarise 
tri_mean

tri_mean$tri_type <- as.factor(tri_mean$tri_type) #factor coercion 

#subset the data: 
type1 <- filter(tri_ct, tri_type == "1")
type2 <- filter(tri_ct, tri_type == "2")


geneA <- expression(paste(bold("wt")))
geneB <- expression(paste(bolditalic("OE-CsHDG5-8")))




p<-ggplot(data = tri_mean, aes(x = tri_type, y=Mcount, fill= tissue)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = TRUE) + 
  geom_errorbar(aes(x = tri_type, ymin = Mcount - SEcount, ymax = Mcount + SEcount, width=0.1),
                position = position_dodge(width=.8), color="black", size=1.0, linetype = 1) +
  scale_fill_discrete(name = "", labels = c(geneA, geneB)) +  
  ylab("Trichome Number") + 
  xlab("") + 
  scale_color_hue(labels = c(geneA, geneB)) + 
  scale_x_discrete(limits = c("1", "2"), 
                   labels = c("Type 1", "Type 2")) + 
  scale_y_continuous(
    limits = c(0,225),
    breaks=seq(0,225, by = 25)) + 
  theme(
    axis.line.y.left = element_line(color = "black"), #color of left axis line
    axis.title.y.left = element_text(color="black", face="bold", size=18), #color, face, and size of y axis title
    axis.text.y.left = element_text(color="black", face="bold", size =18), #color, face, and size of y axis values
    axis.text.x = element_text(color="black", face="bold", size = 18),  #color, face, and size of y axis values
    axis.ticks = element_line(size = 1), #size of tick on x and y axes
    axis.ticks.length = unit(-0.25, "cm"), #makes ticks go inside the plot 
    axis.title.x=element_text(color="black",size=18), #color, size of x axis title (not included)
    plot.title = element_text(face="bold", size=19, hjust = 0.5), 
    plot.caption = element_text(size=12.5, colour ="red"),
    plot.tag = element_text(color="black",face = "bold", size=18),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face ="bold"),
    legend.text = element_text(size = 16, face ="bold"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x =1.2, y = 210, label = "**", size = 8) 
p

fig_7 <- (hdg5.plot|p) #create two-panelled figure 
fig_7

ggsave("06AC_figure.jpeg", device = "jpeg", units="in", width=12, height=5, dpi=300)

