install.packages("scales")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("patchwork")

library("scales")
library(ggplot2)
library(tidyverse)
library(patchwork)

urlfile1="https://raw.githubusercontent.com/apicellap/data/main/CtPharma_GG.csv"
GG <-read.csv(url(urlfile1)) 

urlfile2="https://raw.githubusercontent.com/apicellap/data/main/CW_data_recalc-2.csv"
CW<-read.csv(url(urlfile2)) #data from Cherry Wine (gene expression & metabolite levels)

urlfile3="https://raw.githubusercontent.com/apicellap/data/main/Winter2020_SD_gene.csv"
SD <-read.csv(url(urlfile3)) 

urlfile4="https://raw.githubusercontent.com/apicellap/data/main/Winter2020_wife_gene.csv"
WF<-read.csv(url(urlfile4)) 

GG<-GG %>% rename_with( ~ paste0("GG.", .x)) #adds prefix to every column name in dataframe 
CW<-CW %>% rename_with( ~ paste0("CW.", .x))  
SD<-SD %>% rename_with( ~ paste0("SD.", .x)) 
WF<-WF %>% rename_with( ~ paste0("WF.", .x)) 

GG <- select(GG, 1,18,19 )
names(GG)[names(GG) == 'GG.week'] <- 'week' #changes name of variable from "GG.week" to "week" 
names(GG)[names(GG) == 'GG.HD1.1'] <-  'GG_HDG5' 
names(GG)[names(GG) == 'GG.HD1.se.1'] <- 'GG_HDG5_se' 

GG

CW <- select(CW, 1,26,28 )
names(CW)[names(CW) == 'CW.week'] <- 'week' #changes name of variable
names(CW)[names(CW) == 'CW.HD1.1'] <- 'CW_HDG5' 
names(CW)[names(CW) == 'CW.HD1.se.1'] <- 'CW_HDG5_se' 

CW

SD <- select(SD, 1,12,13 )
names(SD)[names(SD) == 'SD.week'] <- 'week' #changes name of variable
names(SD)[names(SD) == 'SD.HD1.SD'] <- 'SD_HDG5' 
names(SD)[names(SD) == 'SD.HD1.SD.se'] <- 'SD_HDG5_se'
SD

WF <- select(WF, 1,12,13 )
names(WF)[names(WF) == 'WF.week'] <- 'week' #changes name of variable
names(WF)[names(WF) == 'WF.HD1'] <- 'WF_HDG5' 
names(WF)[names(WF) == 'WF.HD1.se'] <- 'WF_HDG5_se' 
WF

vars <- full_join(GG,CW) #combine two dataframes that share a mutual column 
vars <- full_join(vars,SD)
vars <- full_join(vars,WF)
vars


GGp <-ggplot(vars, aes(x = week, y = GG_HDG5)) + geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = week, ymin = GG_HDG5 - GG_HDG5_se, ymax = GG_HDG5 + GG_HDG5_se, width=0.1), size = 1) +
  xlab("") + ylab("") + ggtitle("GG") +
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    limits = c(0,12),
    breaks=c(0,2,4,6,8,10,12)) + 
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
    axis.line = element_line(colour = "black"),
    plot.tag.position = c(0.175,0.85),
    plot.tag = element_text(color="black", face="bold", size = 14, family = "Arial"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x = 2, y = 3.56, label = "*", size = 8, color = "black") +
  annotate("text", x = 3, y = 3.25, label = "*", size = 8, color = "black") +
  annotate("text", x = 4, y = 6.66, label = "*", size = 8, color = "black") +
  annotate("text", x = 5, y = 11.46, label = "*", size = 8, color = "black") +
  annotate("text", x = 6, y =6.38, label = "**", size = 8, color = "black") 
GGp

CWp <-ggplot(vars, aes(x = week, y = CW_HDG5)) + geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = week, ymin = CW_HDG5 - CW_HDG5_se, ymax = CW_HDG5 + CW_HDG5_se, width=0.1), size = 1) +
  xlab("") + ylab("Relative expression") + ggtitle("CW") + 
  scale_x_discrete(name ="", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    limits = c(0,6),
    breaks=c(0,1,2,3,4,5,6)) + 
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
    axis.line = element_line(colour = "black"),
    plot.tag.position = c(0.15,0.85),
    plot.tag = element_text(color="black", face="bold", size = 14, family = "Arial"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x = 5, y = 5.5, label = "**", size = 8, color = "black") 
  
CWp

prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}


SDp <-ggplot(vars, aes(x = week, y = SD_HDG5)) + geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = week, ymin = SD_HDG5 - SD_HDG5_se, ymax = SD_HDG5 + SD_HDG5_se, width=0.1), size = 1) +
  xlab("Week") + ylab("Relative expression") + ggtitle("SD") +  
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    labels = prettyZero, 
    limits = c(0,3),
    breaks=seq(0,3, by = 0.5)) + 
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
    axis.line = element_line(colour = "black"),
    plot.tag.position = c(0.15,0.85),
    plot.tag = element_text(color="black", face="bold", size = 14, family = "Arial"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x = 3, y = 2.8, label = "*", size = 8, color = "black") 
SDp

WFp <-ggplot(vars, aes(x = week, y = WF_HDG5)) + geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = week, ymin = WF_HDG5 - WF_HDG5_se, ymax = WF_HDG5 + WF_HDG5_se, width=0.1), size = 1) +
  xlab("") + ylab("") + ggtitle("WF") +
  scale_x_discrete(name ="", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    limits = c(0,2),
    breaks=seq(0,2, by = 0.25)) + 
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
    axis.line = element_line(colour = "black"),
    plot.tag.position = c(0.15,0.85),
    plot.tag = element_text(color="black", face="bold", size = 14, family = "Arial"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  annotate("text", x = 3, y = 1.025, label = "*", size = 8, color = "black") +
  annotate("text", x = 5, y = 1.9, label = "*", size = 8, color = "black") +
annotate("text", x = 6, y = 0.75, label = "*", size = 8, color = "black") 
  
WFp

pw <- ((CWp|WFp)/(SDp|GGp))
pw

ggsave("03_figure.jpeg", device = "jpeg", units="in", width=10, height=8, dpi=300)

