library(ggplotify)
library(cowplot)
library(ggplot2)
library(patchwork)
library(tidyverse)

#download timecourse data for methyl jasmonate experiment 
urlfile3a="https://raw.githubusercontent.com/apicellap/data/main/autopilot_1h_data.csv"
timecs<-read.csv(url(urlfile3a)) 
head(timecs)

timecs$time <- as.factor(timecs$time) #coerce variables into factors 
timecs$conc <- as.factor(timecs$conc)


### variable - 'conc' ###
# The dataframe was originally written to include 
# four codes for the different concentrations 
# of methyl jasmonate used: 
# conc_a = 0 µM (control)
# conc_b = 100 µM
# conc_c = 500 µM
# conc_d = 1000 µM
# In retrospect, this was unncessary, but 
# I am not going to bother altering 
# the dataframe itself 

#subset data  
tify_100 <- subset(timecs, id2 == "TIFY_AvB") 
tify_500 <- subset(timecs, id2 == "TIFY_AvC")
tify_1000 <- subset(timecs, id2 == "TIFY_AvD")
hdg5_100 <- subset(timecs, id2 == "HD1_AvB")
hdg5_500 <- subset(timecs, id2 == "HD1_AvC")
hdg5_1000 <- subset(timecs, id2 == "HD1_AvD")


# the sole purpose of this graph is to produce the legend
# and convert the legend to an object that is compatible with ggplot2
# The code for this graph also serves as a reference for the rest of the script
# I annotate the lines of code in zz, but not in the rest of the script 
zz <- ggplot(data = timecs, aes(x = time, y=exp, 
                                fill = conc)) + #tells R to fill bar with color according to conc level
  geom_col(width = 0.8, #width of bar
           position = position_dodge2(1), #preserves the horizontal position of bars
           show.legend = TRUE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1), #error bar proportions
                position = position_dodge(width=.8, preserve = "single"), 
                color="black", #error bar color
                size=1.0, #width of error bar
                linetype = 1) + #type of line (solid)
  labs(x="Hours", y="Relative expression") + 
  ggtitle("TIFY9") + 
  scale_fill_manual(values = c("#6abf2d", "#BF392D", "#822DBF", "#2DB3BF"), #pre-selected hexcodes for colors 
                    name= "MeJa\nConcentration (µM)", # \n causes a the following word/phrase to go on the next line
                    labels = c("0", "100", "500","1000")) + #the labels correspond to the hexcode colors
  scale_x_discrete(limits = c("1","3","12" )) + #specify x axis ticks 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black", face="bold",size = 14, family = "Arial"),
    legend.title = element_text(size = 14, face ="bold"),
    legend.text =element_text(size = 12, face ="bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "bold", size = 15)) 
zz

leg <- get_legend(zz)  
leg<-as.ggplot(leg, vjust = 0.4)

geneA1 <- expression(paste(bolditalic("CsTIFY9")))
geneA2 <- expression(paste(bolditalic("CsHDG5")))

a <-ggplot(data = tify_100, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), color="black", size=1.0, linetype = 1) +
  labs(x="", y="Relative expression") + 
  ggtitle(geneA1) + 
  scale_y_continuous(limits = c(0,8),
                     breaks = seq(0,8, by = 1)) +
  scale_fill_manual(values = c("#6abf2d", "#BF392D"), name= "MeJa\nConcentration (µM)", 
                    labels = c("0", "100")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold",size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold"),
    legend.text =element_text(size = 12, face ="bold")) +
  annotate("text", x =0.785, y = 1.5, label = "", size = 7) +
  annotate("text", x = 1.2, y = 6, label = "", size = 7) +
  annotate("text", x = 1.8, y = 3.25, label = "", size = 7) +
  annotate("text", x = 2.20, y = 7.75, label = "", size = 7) +
  annotate("text", x = 2.8, y = 3, label = "", size = 7) +
  annotate("text", x = 3.2, y = 5, label = "", size = 7) 
a

b <-ggplot(data = tify_500, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), color="black", size=1.0, linetype = 1) +
  labs(x="", y="Relative expression") + 
  ggtitle("") + 
  scale_y_continuous(limits = c(0,16),
                     breaks = seq(0,16, by = 4)) +
  scale_fill_manual(values = c("#6abf2d", "#822DBF"), name= "MeJa\nConcentration (µM)", 
                    labels = c("0", "500")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.text =element_text(size = 12, face ="bold"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black",face="bold", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold")) +
  annotate("text", x =0.785, y = 1.5, label = "", size = 7) +
  annotate("text", x = 1.2, y = 10, label = "", size = 7) +
  annotate("text", x = 1.8, y = 3.25, label = "", size = 7) +
  annotate("text", x = 2.215, y = 7.5, label = "", size = 7) +
  annotate("text", x = 2.8, y = 6, label = "", size = 7) +
  annotate("text", x = 3.2, y = 16, label = "", size = 7) 
b



c <-ggplot(data = tify_1000, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), 
                color="black", size=1.0, linetype = 1) +
  labs(x="Hours after MeJA Application", y="Relative expression") + 
  ggtitle("") + 
  scale_y_continuous(limits = c(0,9),
                     breaks = seq(0,9, by = 3)) +
  scale_fill_manual(values = c("#6abf2d","#2DB3BF"),
                    name= "MeJa\nConcentration (µM)", 
                    labels = c("0", "1000")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.text =element_text(size = 12, face ="bold"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black",face="bold", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",face = "bold",size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold")) +
  annotate("text", x =0.785, y = 1.5, label = "", size = 7) +
  annotate("text", x = 1.2, y = 8.5, label = "", size = 7) +
  annotate("text", x = 1.8, y = 3.25, label = "", size = 7) +
  annotate("text", x = 2.215, y = 6.5, label = "", size = 7) +
  annotate("text", x = 2.8, y = 3, label = "", size = 7) +
  annotate("text", x = 3.2, y = 5.5, label = "", size = 7) 
c

d <-ggplot(data = hdg5_100, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), color="black", size=1.0, linetype = 1) +
  labs(x="", y="") + 
  ggtitle(geneA2) + 
  scale_fill_manual(values = c("#6abf2d", "#BF392D"), name= "Concentration \n MeJa (µM)", 
                    labels = c("0", "100")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.text =element_text(size = 12, face ="bold"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold"))
d

e <-ggplot(data = hdg5_500, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), color="black", size=1.0, linetype = 1) +
  labs(x="", y="") + 
  ggtitle("") + 
  scale_y_continuous(limits = c(0,3),
                     breaks = seq(0,3, by = 0.5)) +
  scale_fill_manual(values = c("#6abf2d", "#822DBF"), 
                    name= "Concentration \n MeJa (µM)", 
                    labels = c("0","500")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.text =element_text(size = 12, face ="bold"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black",face="bold", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",size = 14, face="bold",family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold"))
e

f <-ggplot(data = hdg5_1000, aes(x = time, y=exp, fill = conc)) +
  geom_col(width = 0.8, position = position_dodge2(1), show.legend = FALSE) + 
  geom_errorbar(aes(x=time, ymin=exp *1 -exp.se *1, ymax=exp*1+exp.se, width=.1),
                position = position_dodge(width=.8, preserve = "single"), color="black", size=1.0, linetype = 1) +
  labs(x="Hours after MeJA Application", y="") + 
  ggtitle("") + 
  scale_y_continuous(limits = c(-0.1,2),
                     breaks = seq(0,2, by = 0.5)) +
  scale_fill_manual(values = c("#6abf2d", "#2DB3BF"), name= "Concentration \n MeJa (µM)", 
                    labels = c("0", "1000")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.text =element_text(size = 12, face ="bold"),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black",face="bold", size = 14, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.text.x = element_text(color="black", face="bold", size = 14, family = "Arial"),
    axis.title.x=element_text(color="black",face="bold",size = 14, family = "Arial"),
    plot.title = element_text(hjust = 0.5, family = "Arial", face = "italic", size = 15),
    legend.title = element_text(size = 14, face ="bold"))
f


PW_1h <- (((a|d)/(b|e)/(c|f)))/leg

PW_1h + plot_layout(guides = "collect") &
  theme(plot.title = element_text(hjust = 0.5, family = "Arial", face="bold", size=14),
        plot.tag = element_text(family = "Arial", face = "bold", size = 14),
        legend.text = (element_text(size = 12, face = "bold"))) 


# the final plot may be cropped to remove white space at bottom of figure:
ggsave("APtimecs.jpeg", device = "jpeg", units="in", width=8.5, height=12, dpi=175) 
