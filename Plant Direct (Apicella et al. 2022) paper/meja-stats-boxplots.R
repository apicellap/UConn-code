#install/load packages: 
install.packages("ggplot2")
install.packages("agricolae")
install.packages("dplyr")
install.packages("patchwork")
install.packages("ggpubr")
library(ggplot2)
library(agricolae) 
library(tidyverse)
library(patchwork)
library(ggpubr)

# download the data: 
urlfile3="https://raw.githubusercontent.com/apicellap/data/main/CtPharma_dwt.csv"
MJ<-read.csv(url(urlfile3))
head(MJ)

# Dataframe modifications: 
names(MJ)[names(MJ) == 'X.DWT'] <- 'perc.dwt' #changes name of variable
MJ$dec.dwt <- (MJ$perc.dwt/100) #creates new column where percentage is converted to decimal, 
                                #required by arcsin() calculation
MJ$ASin <-asin(sqrt(MJ$dec.dwt)) #arcsin calculation. 
MJ$CONC <- as.factor(MJ$CONC) #making these variables factors
MJ$WEEK <- as.factor(MJ$WEEK)
MJ$METAB <- as.factor(MJ$METAB)
MJ <- mutate(MJ, CONC = relevel(CONC, ref = "0")) #make 0 µM the reference level 

### subset data: 
THC <- subset(MJ, METAB=="THC")
TOTAL <- subset(MJ, METAB=="TOTAL")
THC_wk1 <- subset(THC, WEEK=="1")
THC_wk2 <- subset(THC, WEEK=="2")
THC_wk3 <- subset(THC, WEEK=="3")
THC_wk4 <- subset(THC, WEEK=="4")
TOTAL_wk1 <- subset(TOTAL, WEEK=="1")
TOTAL_wk2 <- subset(TOTAL, WEEK=="2")
TOTAL_wk3 <- subset(TOTAL, WEEK=="3")
TOTAL_wk4 <- subset(TOTAL, WEEK=="4")

### create linear models, perform ANOVA, and Fisher's LSD test: 
### Linear models for THC levels in different weeks
mod1 <- lm(ASin ~ CONC, THC_wk1) #linear model for the dataframe THC_wk1, does ASin change as a function of CONC?
summary(mod1) #ANOVA
LSD1 <- LSD.test(mod1, "CONC", group=TRUE) #perform LSD test
LSD1 #print LSD results 

#the next six lines help me automate adding the LSD letters to one plot but this isn't necessary 
CONC <- as.list(c(1000,0,500,100))
df_new<-LSD1$groups 
df_new$CONC <- CONC
df_new$CONC <- as.numeric(df_new$CONC)
df_new <- arrange(df_new, CONC)
df_new %>% mutate_if(is.character, str_to_upper) -> df_new

mod2 <- lm(ASin ~ CONC, THC_wk2)
summary(mod2)
LSD2 <- LSD.test(mod2, "CONC", group=TRUE)
LSD2
P2<-plot(LSD2, xlab="Concentration", ylab="THC", main="week 2") 

mod3 <- lm(ASin ~ CONC, THC_wk3)
summary(mod3)
LSD3 <- LSD.test(mod3, "CONC", group=TRUE)
LSD3
P3<-plot(LSD3, xlab="Concentration", ylab="THC", main="week 3") 

mod4 <- lm(ASin ~ CONC, THC_wk4)
summary(mod4)
LSD4 <- LSD.test(mod4, "CONC", group=TRUE)
LSD4
P4<-plot(LSD4, xlab="Concentration", ylab="THC", main="week 4") 

### linear models for total cannabinoid analysis 

mod5 <- lm(ASin ~ CONC, TOTAL_wk1)
summary(mod5)
LSD5 <- LSD.test(mod5, "CONC", group=TRUE)
LSD5
plot(LSD5, xlab="Concentration", ylab="TOTAL", main="week 1") 

mod6 <- lm(ASin ~ CONC, TOTAL_wk2)
summary(mod6)
LSD6 <- LSD.test(mod6, "CONC", group=TRUE)
LSD6
plot(LSD6, xlab="Concentration", ylab="TOTAL", main="week 2") 

mod7 <- lm(ASin ~ CONC, TOTAL_wk3)
summary(mod7)
LSD7 <- LSD.test(mod7, "CONC", group=TRUE)
LSD7
plot(LSD7, xlab="Concentration", ylab="TOTAL", main="week 3") 

mod8 <- lm(ASin ~ CONC, TOTAL_wk4)
summary(mod8)
LSD8 <- LSD.test(mod8, "CONC", group=TRUE)
LSD8
plot(LSD8, xlab="Concentration", ylab="TOTAL", main="week 4") 


###### Boxplots #####

colors <- c("#00AFBB", "#E7B800", "#808080", "#FC4E07") #hexcodes for colors of boxplots
names(colors) = c("0", "100", "500", "1000") #gives names to the colors 

w1THC <-ggboxplot(THC_wk1, x = "CONC", y = "perc.dwt", 
                  color = "CONC", 
                  palette = colors, 
                  add = "jitter", 
                  xlab = "" , 
                  ylab = "THC content (% weight)",
                  title = "Week 1",
                  ylim = c(2,3.5),
                  yticks.by = 0.25,
                  legend.title = "Concentration (µM)") +  rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.y.left = element_line(color = "black"),
        axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
        axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial")) + 
  annotate("text", x = 1:4, y = c(3.25,2.8,3,3.2), 
           label = df_new$groups, #adding labels like this takes work, but if the analysis changes, this can remain the same and will automatically update
           size = 6, fontface = "bold")
w1THC


w2THC <-ggboxplot(THC_wk2, x = "CONC", y = "perc.dwt", 
                  color = "CONC", 
                  palette = colors, 
                  add = "jitter", 
                  ylim = c(3,6),
                  xlab = "" , 
                  ylab = "",
                  title = "Week 2") + rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    # axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 5.25, label = "C", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 5.25, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 5, label = "BC", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 4, y = 5.75, label = "AB", size = 6, fontface = "bold", family= theme_get()$text[["family"]])

w2THC

# ggsave("meja_THC-wk2.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


w3THC <-ggboxplot(THC_wk3, x = "CONC", y = "perc.dwt", 
                  color = "CONC", 
                  palette = colors, 
                  add = "jitter", 
                  ylim = c(6,11),
                  xlab = "MeJA Concentration (µM)" , 
                  ylab = "THC Content (% weight)",
                  title = "Week 3") +  rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 10.5, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 10, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 10.5, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) 

w3THC

# ggsave("meja_THC-wk3.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)

w4THC <-ggboxplot(THC_wk4, x = "CONC", y = "perc.dwt", 
                  color = "CONC", 
                  palette = colors, 
                  add = "jitter", 
                  ylim = c(6,8),
                  xlab = "MeJA Concentration (µM)" , 
                  ylab = "",
                  title = "Week 4") + rremove("legend") + 
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 7.5, label = "B", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 7.75, label = "AB", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 7.9, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) 

w4THC
# ggsave("meja_THC-wk4.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


THC_meja_PW <- ((w1THC|w2THC)/(w3THC|w4THC))
THC_meja_PW 

ggsave("thc-meja_PW.jpeg", device = "jpeg", units="in", width=10, height=8, dpi=300)


w1TOTAL <-ggboxplot(TOTAL_wk1, x = "CONC", y = "perc.dwt", 
                    color = "CONC", 
                    palette = colors, 
                    add = "jitter", 
                    xlab = "" , 
                    ylab = "Total Cannabinoid Content \n (% weight)",
                    title = "Week 1",
                    ylim = c(2.5,5),
                    yticks.by = 0.25,
                    legend.title = "Concentration (µM)") +  rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 4.6, label = "B", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 4.1, label = "B", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 4.25, label = "B", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 4, y = 4.6, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]])


w1TOTAL

# ggsave("meja_TOTAL-wk1.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


w2TOTAL <-ggboxplot(TOTAL_wk2, x = "CONC", y = "perc.dwt", 
                    color = "CONC", 
                    palette = colors, 
                    add = "jitter", 
                    ylim = c(4,8),
                    xlab = "" , 
                    ylab = "",
                    title = "Week 2") + rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 7, label = "C", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 7.25, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y =6.75, label = "BC", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 4, y = 7.6, label = "AB", size = 6, fontface = "bold", family= theme_get()$text[["family"]])

w2TOTAL

# ggsave("meja_TOTAL-wk2.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


w3TOTAL <-ggboxplot(TOTAL_wk3, x = "CONC", y = "perc.dwt", 
                    color = "CONC", 
                    palette = colors, 
                    add = "jitter", 
                    ylim = c(9,13),
                    xlab = "MeJA Concentration (µM)" , 
                    ylab = "Total Cannabinoid Content \n (% weight)",
                    title = "Week 3") +  rremove("legend") +
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 12.75, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 11.5, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 11, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) 

w3TOTAL

# ggsave("meja_TOTAL-wk3.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


w4TOTAL <-ggboxplot(TOTAL_wk4, x = "CONC", y = "perc.dwt", 
                    color = "CONC", 
                    palette = colors, 
                    add = "jitter", 
                    ylim = c(8,10.5),
                    xlab = "MeJA Concentration (µM)" , 
                    ylab = "",
                    title = "Week 4") + rremove("legend") + 
  theme(plot.title = element_text(size= 14, hjust = 0.5, family = "Arial", face = "bold")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.title.y.left = element_text(color="black", face="bold", size = 16, family = "Arial"),
    axis.text.y.left = element_text(color="black", face="bold", size = 16, family = "Arial",margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.x = element_text(color="black", face="bold", size = 16, family = "Arial", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x=element_text(color="black", face="bold",size = 16, family = "Arial"))+ 
  annotate("text", x = 1, y = 10, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 2, y = 10, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) +
  annotate("text", x = 3, y = 10, label = "A", size = 6, fontface = "bold", family= theme_get()$text[["family"]]) 

w4TOTAL
# ggsave("meja_TOTAL-wk4.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)


TOTAL_meja_PW <- ((w1TOTAL|w2TOTAL)/(w3TOTAL|w4TOTAL))
TOTAL_meja_PW
ggsave("total-meja.jpeg", device = "jpeg", units="in", width=10, height=8, dpi=300)

