install.packages("patchwork")
library(patchwork)
library(tidyverse)
library(stringr)
library(ggplot2)

urlfile1f="https://raw.githubusercontent.com/apicellap/data/main/fig8a.csv"
ge<-read.csv(url(urlfile1f)) 

geneA1 <- expression(paste(bold("Control")))
gene1 <- expression(paste(bolditalic("CsHDG5")))


ge_plot <-ggplot(ge, aes(x = treat, y = exp)) +
  geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = treat, ymin = exp - sd_m, ymax = exp + sd_p, width=0.1), size = 1) +
  xlab("") + ylab("Relative expression")  + ggtitle(gene1)+
  scale_x_discrete(limits = c("ctl", "atpep3"), 
                   labels = c("Control", "AtPep3")) + 
  scale_y_continuous(
     limits = c(0,50),
    breaks= seq(0,50, by = 10)) + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
   annotate("text", x = 2, y = 47, label = "***", size = 8, color = "black") 
ge_plot

###### cleaning up data that is not formatted to R's liking 

urlfile1g="https://raw.githubusercontent.com/apicellap/data/main/nasty_ros.csv"
nas_ros<-read.csv(url(urlfile1g)) 
View(nas_ros)

#data wrangling: (the following code must be run in order and each line must be ran exactly once)
#start from the first line in this series if you get an error message
#in the console, you may want to periodically check what the modifications do to the dataframes
  #this can be accomplished by entering the head(<insert name of dataframe>) and 
  #View(<insert name of dataframe>) commands in the console 
wt <- select(nas_ros, 1:7,20) #isolate columns associated with wild type data in one dataframe 
wt <- gather(wt, arb, ros,2:7,na.rm = FALSE, convert = FALSE) 
wt <- select(wt, X, ros,sec)
names(wt)[names(wt) == 'X'] <- 'cycle' #changes name of variable
names(wt)[names(wt) == 'wt'] <- 'ros' #changes name of variable
wt$treat <- "wt" #create new variable in dataframe and add the value in quotes to every observation

trt <- select(nas_ros, 1,11:16,20) #select only certain numbered columns of the original dataframe
trt <- gather(trt, arb, ros,2:7,na.rm = FALSE, convert = FALSE)
head(trt)
names(trt)[names(trt) == 'X'] <- 'cycle' #changes name of variable
names(trt)[names(trt) == 'trt'] <- 'ros' #changes name of variable
trt <- select(trt, cycle, ros, sec)
trt$treat <- "trt"

tidyros<-rbind(wt, trt)
tidyros$sec <- as.numeric(tidyros$sec)
tidyros <- filter(tidyros, sec <= 3424)
tidyros$min <- tidyros$sec/60
tidyros %>% mutate(across(starts_with("min"), round, 1))
tidyros$min <- as.factor(tidyros$min)

sum <- tidyros %>% 
  group_by(cycle, treat, min, .groups = TRUE) %>%
  summarise(
    meanT = mean(ros),
    se = (sd(ros)/sqrt(6))) 
head(sum)

sum <- select(sum, 1:3,5,6)
head(sum)

sum$min <- as.numeric(sum$min)
sum %>% mutate(across(starts_with("min"), round, 2))


sum <- arrange(sum, by_group = treat)
View(sum)


sum$treat <- gsub("trt","OE-CsHDG5-8" , sum$treat)
sum$treat <- as.factor(sum$treat)

sum <- mutate(sum,treat = relevel(treat, ref = "wt")) # the wt level had to be coerced into the reference level 
#this allows me to assign it to a different line color 
levels(sum$treat)

sum$meanT <- sum$meanT/100
sum$se   <- sum$se/100
View(sum)

#############################################################################################

geneA <- expression(paste(bold("wt")))
geneB <- expression(paste(bolditalic("OE-CsHDG5-8")))

rosline <- ggplot(sum, aes(x = min, y = meanT, color = treat)) + 
  geom_errorbar(aes(x = min, ymin = meanT - se, ymax = meanT + se)) +
  geom_point(size =4) + 
  geom_line(mapping = aes(x=min, y = meanT, color = treat), size =2) +
  ylab("RLUs (x100)") + xlab("Time (minutes)") + ggtitle("AtPep3") + 
  scale_fill_continuous(name = "", labels = c(geneA, geneB)) +  
  scale_x_continuous(limits = c(1,33),
                     breaks = seq(1, 33, by = 2)) +
scale_y_continuous(limits = c(0,28),
                      breaks = seq(0,28, by = 4)) +
  scale_color_hue(labels = c(geneA, geneB)) + 
theme(
    axis.line.y.left = element_line(color = "black"), #color of left axis line
    axis.title.y.left = element_text(color="black", face="bold", size=18), #color, face, and size of y axis title
    axis.text.y.left = element_text(color="black", face="bold", size =18), #color, face, and size of y axis values
    axis.text.x = element_text(color="black", face="bold", size =14, vjust=0.5),  #color, face, and size of y axis values
    axis.ticks = element_line(size = 1), #size of tick on x and y axes
    axis.ticks.length = unit(-0.25, "cm"), #makes ticks go inside the plot 
    axis.title.x=element_text(color="black",face = "bold", size=18, vjust =0.6), #color, size of x axis title (not included)
    plot.title = element_text(face="bold", size=19, hjust = 0.5), 
    plot.caption = element_text(size=12.5, colour ="red"),
    plot.tag = element_text(color="black",face = "bold", size=18),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face ="bold"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) 

rosline

fig8 <- (ge_plot|rosline)
fig8
fig8 + 
  plot_layout(widths = c(1, 2.5)) #this makes the right panel have a greater proportion of the figure than the left panel 

ggsave("08_figure.jpeg", device = "jpeg", units="in", width=12, height=5, dpi=300)



