#install/load packages: 
installed.packages("ggplot2")
installed.packages("tidyverse")
installed.packages("patchwork")

library(ggplot2)
library(tidyverse)
library(patchwork)

#download data: 
urlfile1b="https://raw.githubusercontent.com/apicellap/data/main/CtPharma_GG.csv"
gg_data<-read.csv(url(urlfile1b)) #data from CtPharma Gorilla Glue (gene expression & metabolite levels)
urlfile2b="https://raw.githubusercontent.com/apicellap/data/main/CW_data_recalc-2.csv"
cw_data<-read.csv(url(urlfile2b)) #data from Cherry Wine (gene expression & metabolite levels)


#dataframe modifications: 
gg_data<-gg_data %>% rename_with( ~ paste0("GG.", .x)) #adds prefix to every column name in dataframe 
cw_data<-cw_data %>% rename_with( ~ paste0("CW.", .x)) 
names(gg_data)[names(gg_data) == 'GG.week'] <- 'week' #changes name of variable
names(cw_data)[names(cw_data) == 'CW.week'] <- 'week'
GG.CW=full_join(gg_data, cw_data) #better way to combine two dataframes. At least one variable must be in common 
                                                  #these two dfs were combined since they share the 'week' variable 
###################################################################################

#Plot the data for PT1/4 genes with CBG content: 
Prim.y_title.001 <- expression(paste("Relative expression of ", italic("PT4")))
ggtitle.001 <- expression(paste(italic("PT4"), " expression vs. CBG content in GG"))

GG.CBG.PT4<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.PT4 *1), stat="identity",fill="#ba5b2f", color = "#6b3d25", size = 1) +
  geom_errorbar(aes(x=week, ymin=GG.PT4 *1 -GG.PT4.se *1, ymax=GG.PT4*1, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=GG.CBG*(4), group=1), size=1) +
  geom_errorbar(aes(x=week, ymin=GG.CBG *4, ymax=GG.CBG *4 +GG.CBG.se *4, width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  geom_point(size = 2.75,aes(x=week, y=GG.CBG*(4), group=1)) +
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,2,4,6,8), limits =c(-1,8),
    name= Prim.y_title.001,
    sec.axis = sec_axis(~ . *(1/4), name="CBG content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(face="bold", size=19,hjust = 0.5)) +
   ggtitle(ggtitle.001) +
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 4, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 6.5, label = "*", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = .5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = .5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y =.51, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 1, label = "", size = 8, color = "#6b3d25")
GG.CBG.PT4

Prim.y_title.002 <- expression(paste("Relative expression of ", italic("PT1")))
ggtitle.002 <- expression(paste(italic("PT1"), " expression vs. CBG content in GG"))


GG.CBG.PT1<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.PT1 *1), stat="identity",fill="#ba5b2f", color = "#6b3d25", size = 1) +
  geom_errorbar(aes(x=week, ymin=GG.PT1 *1 -GG.PT1.se *1, ymax=GG.PT1*1, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=GG.CBG*(21), group=1), size=1) +
  geom_errorbar(aes(x=week, ymin=GG.CBG *21, ymax=GG.CBG *21 +GG.CBG.se *21, width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  geom_point(size = 2.75,aes(x=week, y=GG.CBG*(21), group=1)) +
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,10,20,30,40), limits =c(-1,42),
    name= Prim.y_title.002,
    sec.axis = sec_axis(~ . *(1/21), name="CBG content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.002) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 2, y = 1, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 15, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = 39, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = 22, label = "*", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y =1, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 1, label = "", size = 8, color = "#6b3d25")
GG.CBG.PT1


Prim.y_title.003 <- expression(paste("Relative expression of ", italic("PT4")))
ggtitle.003 <- expression(paste(italic("PT4"), " expression vs. CBG content in CW"))


CW.PT4.CBG<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.PT4 *1), stat="identity",fill="#2f95ba", color = "#25516b", size = 1) +
  geom_errorbar(aes(x=week, ymin=CW.PT4 *1 -CW.PT4.se *1, ymax=CW.PT4*1, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=CW.CBG*(65), group=1), size=1) +
  geom_errorbar(aes(x=week, ymin=CW.CBG *65, ymax=CW.CBG *65 +CW.CBG.se *65, width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  geom_point(size = 2.75,aes(x=week, y=CW.CBG*(65), group=1)) +
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(-2,0,3,6,9,12,15), limits =c(-2,10),
    name= Prim.y_title.003,
    sec.axis = sec_axis(~ . *(1/65), name="CBG content (% dry weight)"))+
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18)) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.003) +
  theme(plot.title = element_text(face="bold", size=19)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 1.25, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 4, y = 9, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 5, y = 3, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 6, y = 2.5, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 7.5, label = "", size = 8, color = "#25516b")

CW.PT4.CBG


Prim.y_title.004 <- expression(paste("Relative expression of ", italic("PT1")))
ggtitle.004 <- expression(paste(italic("PT1"), " expression vs. CBG content in CW"), col.main = "black")


CW.PT1.CBG<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.PT1 *1), stat="identity",fill="#2f95ba", color = "#25516b", size = 1) +
  geom_errorbar(aes(x=week, ymin=CW.PT1 *1 -CW.PT1.se *1, ymax=CW.PT1*1, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  #creates linegraph/connects points from geom_point
  geom_line(aes(x=week, y=CW.CBG*(60), group=1), size =1) +
  geom_errorbar(aes(x=week, ymin=CW.CBG *60, ymax=CW.CBG *60 +CW.CBG.se *60, width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  geom_point(size = 2.75, aes(x=week, y=CW.CBG*(60), group=1)) +
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,1,2,3,4,5,6,7,8,9), limits =c(0,9),
    name= Prim.y_title.004,
    sec.axis = sec_axis(~ . *(1/60), name="CBG content (% dry weight)"))+
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.004) +
  theme(plot.title = element_text(face="bold", size=19,color="black")) +
  theme(plot.title = element_text(hjust = 0.5,)) + 
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 1.25, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8,color = "#25516b") +
  annotate("text", x = 4, y = 3, label = "", size = 8,color = "#25516b") +
  annotate("text", x = 5, y = 8, label = "*", size = 8,color = "#25516b") +
  annotate("text", x = 6, y = 3.5, label = "*", size = 8,color = "#25516b") +
  annotate("text", x = 7, y = 6.5, label = "", size = 8,color = "#25516b")

CW.PT1.CBG

GG.CW.PW <- (CW.PT1.CBG/GG.CBG.PT1)|(CW.PT4.CBG/GG.CBG.PT4)
GG.CW.PW 
ggsave("01_figure.jpeg", device = "jpeg", units="in", width=10, height=8, dpi=300)

###################################################################################

#Plot the data for CBDAS expression/CBD content & THCAS expression/THC content
Prim.y_title.005 <- expression(paste("Relative expression of ", italic("CBDAS")))
ggtitle.005 <- expression(paste(italic("CBDAS"), " expression vs. CBD content in CW"))

CW.CBDAS.CBD<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.CBDAS *1), stat="identity",fill="#2f95ba", color = "#25516b", size = 1) +
  geom_errorbar(aes(x=week, ymin=CW.CBDAS *1 -CW.CBDAS.se *1, ymax=CW.CBDAS*1, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=CW.CBD*(1/3), group=1), size=1) +
  geom_point(size = 2.75,aes(x=week, y=CW.CBD*(1/3), group=1)) +
  geom_errorbar(aes(x=week, ymin=CW.CBD*(1/3), ymax=CW.CBD*(1/3)+CW.CBD.se*(1/3), width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,2,4,6,8), limits =c(0,8),
    name= Prim.y_title.005,
    sec.axis = sec_axis(~ . *(3), name="CBD content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.005) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 2, y = 1, label = "*", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 4, y = 2.4, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 5, y = 6, label = "*", size = 8, color = "#25516b") +
  annotate("text", x = 6, y = 4.5, label = "*", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 7.5, label = "", size = 8, color = "#25516b")

CW.CBDAS.CBD


Prim.y_title.006 <- expression(paste("Relative expression of ", italic("THCAS")))
ggtitle.006 <- expression(paste(italic("THCAS"), " expression vs. THC content in CW"))

CW.THCAS.THC<-ggplot(data=GG.CW
                     ) +
  geom_bar(mapping= aes(x=week, y=CW.THCAS *1), stat="identity",fill="#2f95ba", color = "#25516b", size = 1) +
  geom_errorbar(aes(x=week, ymin=CW.THCAS *1 -CW.THCAS.se *1, ymax=CW.THCAS*1, width=.2),
                position = position_dodge(0.9), color="#25516b", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=CW.THC*10, group=1), size=1) +  #creates linegraph/connects points from geom_point
  geom_point(size = 2.75, aes(x=week, y=CW.THC*10, group=1)) +
  geom_errorbar(aes(x=week, ymin=CW.THC *10, ymax=CW.THC*10+CW.THC.se*10, width=.2),
                position = position_dodge(0.9), color="black", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,2,4,6,8,10), limits =c(-1,10),
    name= Prim.y_title.006,
    sec.axis = sec_axis(~ . *0.1, name="THC content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color="#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color="#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18)) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.006) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 2, y = 2.5, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 1, label = "", size = 8,color = "#25516b") +
  annotate("text", x = 4, y = 2.5, label = "", size = 8,color = "#25516b") +
  annotate("text", x = 5, y = 6, label = "", size = 8,color = "#25516b") +
  annotate("text", x = 6, y = 8, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 8.5, label = "", size = 8,color = "#25516b")

CW.THCAS.THC

Prim.y_title.007 <- expression(paste("Relative expression of ", italic("CBDAS")))
ggtitle.007 <- expression(paste(italic("CBDAS"), " expression vs. CBD content in GG"))

GG.CBDAS.CBD<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.CBDAS *1), stat="identity",fill="#ba5b2f", color = "#6b3d25", size = 1) +
  geom_errorbar(aes(x=week, ymin=GG.CBDAS *1 -GG.CBDAS.se *1, ymax=GG.CBDAS*1, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=GG.CBD*(10), group=1), size=1) +
  geom_point(size = 2.75, aes(x=week, y=GG.CBD*(10), group=1)) +
  geom_errorbar(aes(x=week, ymin=GG.CBD *10, ymax=GG.CBD*10+GG.CBD.se*10, width=.2),
                position = position_dodge(0.9), color = "black", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,0.5,1,1.5,2,2.5), limits =c(0,2.5),
    name= Prim.y_title.007,
    sec.axis = sec_axis(~ . *(1/10), name="CBD content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.007) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 2, y = 1.5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 1.5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = 1.5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = 1.5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y = 0.6, label = "*", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 1.9, label = "", size = 8, color = "#6b3d25")
GG.CBDAS.CBD

Prim.y_title.008 <- expression(paste("Relative expression of ", italic("THCAS")))
ggtitle.008 <- expression(paste(italic("THCAS"), " expression vs. THC content in GG"))

GG.THC.THCAS<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.THCAS *1), stat="identity",fill="#ba5b2f", color = "#6b3d25", size = 1) +
  geom_errorbar(aes(x=week, ymin=GG.THCAS *1 -GG.THCAS.se *1, ymax=GG.THCAS*1, width=.2),
                position = position_dodge(0.9), color="#6b3d25", size=1.0, linetype = 1)+
  geom_line(aes(x=week, y=GG.THC*(1/20), group=1), size=1) +
  geom_point(size = 2.75, aes(x=week, y=GG.THC*(1/20), group=1)) +
  geom_errorbar(aes(x=week, ymin=GG.THC *(1/20), ymax=GG.THC*(1/20)+GG.THCAS.se*(1/20), width=.2),
                position = position_dodge(0.9), color="black", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,0.2,0.4,0.6,0.8,1,1.2), limits =c(0,1.2),
    name= Prim.y_title.008,
    sec.axis = sec_axis(~ . *(20), name="THC content (% dry weight)")) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color="#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color="#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18)) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.008) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate("text", x = 2, y = 1.2, label = "", size = 8) +
  annotate("text", x = 3, y = 1.2, label = "", size = 8) +
  annotate("text", x = 4, y = 1.2, label = "", size = 8) +
  annotate("text", x = 5, y = 1.2, label = "", size = 8) +
  annotate("text", x = 6, y = 1.2, label = "", size = 8) +
  annotate("text", x = 7, y = 1.2, label = "", size = 8)
GG.THC.THCAS


GG.CW.F2.pw <- ((CW.CBDAS.CBD/GG.CBDAS.CBD)|(CW.THCAS.THC/GG.THC.THCAS))
GG.CW.F2.pw 
ggsave("02_figure.jpeg", device = "jpeg", units="in", width=12, height=8, dpi=300)

###################################################################################

#Plot the data for GPPS, TKS, and OAC expression levels: 
Prim.y_title.009 <- expression(paste("Relative expression of ", italic("GPPS")))
ggtitle.009 <- expression(paste(italic("GPPS"), " expression in CW"))


CW.GPPS<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.GPPS *1), stat="identity",fill="#2f95ba") +
  geom_errorbar(aes(x=week, ymin=CW.GPPS *1 -CW.GPPS.se *1, ymax=CW.GPPS*1+CW.GPPS.se, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,2,4,6,8,10), limits =c(-0.5,10),
    name= Prim.y_title.009) +
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.009) +
  theme(plot.title = element_text(face="bold", size=19)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 4, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 5, y = 9.5, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 6, y = 6.5, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 1.5, label = "", size = 8, color = "#25516b")

CW.GPPS

Prim.y_title.010 <- expression(paste("Relative expression of ", italic("TKS")))
ggtitle.010 <- expression(paste(italic("TKS"), " expression in CW"))


CW.TKS<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.TKS *1), stat="identity",fill="#2f95ba") +
  geom_errorbar(aes(x=week, ymin=CW.TKS *1 -CW.TKS.se *1, ymax=CW.TKS*1+CW.TKS.se, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,1,2,3,4,5,6,7), limits =c(0,7),
    name= Prim.y_title.010) +
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.010) +
  theme(plot.title = element_text(face="bold", size=19)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 4, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 5, y = 6, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 6, y = 6.75, label = "**", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 1.5, label = "", size = 8, color = "#25516b")

CW.TKS


Prim.y_title.011 <- expression(paste("Relative expression of ", italic("GPPS")))
ggtitle.011 <- expression(paste(italic("GPPS"), " expression in GG"))


GG.GPPS<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.GPPS *1), stat="identity",fill="#ba5b2f") +
  geom_errorbar(aes(x=week, ymin=GG.GPPS *1 -GG.GPPS.se *1, ymax=GG.GPPS*1+GG.GPPS.se, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,3,6,9,12,15,18,21), limits =c(-2,21),
    name= Prim.y_title.011) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18)) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.011) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Week") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2, y = 8, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = 8, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = 8, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y = 21, label = "*", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 8, label = "", size = 8, color = "#6b3d25")
GG.GPPS

Prim.y_title.012 <- expression(paste("Relative expression of ", italic("TKS")))
ggtitle.012 <- expression(paste(italic("TKS"), " expression in GG"))


GG.TKS<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.TKS *1), stat="identity",fill="#ba5b2f") +
  geom_errorbar(aes(x=week, ymin=GG.TKS *1 -GG.TKS.se *1, ymax=GG.TKS*1+GG.TKS.se, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,1,2,3,4,5), limits =c(-0.25,5),
    name= Prim.y_title.012) + 
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18),
        axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_text(color="black",size=18)) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.012) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Weeks after Flower Initiation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y = 1, label = "***", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 2.5, label = "***", size = 8, color = "#6b3d25")
GG.TKS


Prim.y_title.013 <- expression(paste("Relative expression of ", italic("OAC")))
ggtitle.013 <- expression(paste(italic("OAC"), " expression in GG"))


GG.OAC.nb<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=GG.OAC *1), stat="identity",fill="#ba5b2f") +
  geom_errorbar(aes(x=week, ymin=GG.OAC *1 -GG.OAC.se *1, ymax=GG.OAC*1+GG.OAC.se, width=.2),
                position = position_dodge(0.9), color = "#6b3d25", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,50,100,150,200,250,300,350,400), limits =c(-25,410),
    name= Prim.y_title.013) + 
  theme(
    axis.line.y.right = element_line(color = "black"),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color = "#6b3d25", face="bold", size=18, vjust = -1,),
    axis.text.y.left = element_text(color = "#6b3d25", face="bold", size =18),
    axis.text.x = element_text(color="black", face="bold", size = 18),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.ticks.y.right = element_blank(),
    axis.title.x=element_text(color="black",size=18),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.013) +
  theme(plot.title = element_text(face="bold", size=19)) +
  xlab("Week") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 3, y = 5, label = "", size = 8, color = "#6b3d25") +
  annotate("text", x = 4, y = 55, label = "**", size = 8, color = "#6b3d25") +
  annotate("text", x = 5, y = 55, label = "*", size = 8, color = "#6b3d25") +
  annotate("text", x = 6, y = 77, label = "***", size = 8, color = "#6b3d25") +
  annotate("text", x = 7, y = 360, label = "**", size = 8, color = "#6b3d25")  


GG.OAC.nb 


Prim.y_title.014 <- expression(paste("Relative expression of ", italic("OAC")))
ggtitle.014 <- expression(paste(italic("OAC"), " expression in CW"))


CW.OAC<-ggplot(data=GG.CW) +
  geom_bar(mapping= aes(x=week, y=CW.OAC *1), stat="identity",fill="#2f95ba") +
  geom_errorbar(aes(x=week, ymin=CW.OAC *1 -CW.OAC.se *1, ymax=CW.OAC*1+CW.OAC.se, width=.2),
                position = position_dodge(0.9), color = "#25516b", size=1.0, linetype = 1)+
  scale_x_discrete(name ="Week", 
                   limits=c("1","2","3", "4", "5", "6", "7")) +
  scale_y_continuous(
    breaks=c(0,2,4,6,8,10,12,14,16), limits =c(0,16),
    name= Prim.y_title.014) +
  # sec.axis = sec_axis(~ . *(1/100), name="CBG content (% dry weight)"))+
  theme(axis.line.y.right = element_line(color = "black"),
        axis.title.y.right = element_text(color="black", size = 18),
        axis.text.y.right = element_text(color="black", face="bold", size = 18),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color = "#25516b", face="bold", size=18),
        axis.text.y.left = element_text(color = "#25516b", face="bold", size =18),
        axis.text.x = element_text(color="black", face="bold", size = 18),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.title.x=element_blank()) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  ggtitle(ggtitle.014) +
  theme(plot.title = element_text(face="bold", size=19)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Weeks after Flower Initiation") +
  annotate("text", x = 2, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 3, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 4, y = 2, label = "", size = 8, color = "#25516b") +
  annotate("text", x = 5, y = 7, label = "*", size = 8, color = "#25516b") +
  annotate("text", x = 6, y = 7, label = "****", size = 8, color = "#25516b") +
  annotate("text", x = 7, y = 15.75, label = "*", size = 8, color = "#25516b")

CW.OAC

pw <- ((CW.GPPS/GG.GPPS)|(CW.TKS/GG.TKS)|(CW.OAC/GG.OAC.nb))
pw 

ggsave("gpps-oac-tks.jpeg", device = "jpeg", units="in", width=14, height=9, dpi=300)