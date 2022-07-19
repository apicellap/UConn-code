
library(ggplot2)
library(dplyr)
library(patchwork)

urlfile1c="https://raw.githubusercontent.com/apicellap/data/main/hd1_tissuetype_raw.csv"
hd1<-read.csv(url(urlfile1c)) 
head(hd1)

hd1 <- subset(hd1, target_gene == "HD1") #subset dataframe to remove rows with missing values 


# qPCR data analysis: 
ref_df<-filter(hd1, tissue == "root") #subset tissue by root to isolate and do calculations on its mean ct values 
Gmean_ct_target <- mean(ref_df$target_mean_ct) #calculates means of the mean ct values in the target gene 
Gmean_ct_hk <- mean(ref_df$hk_mean_ct)   #calculates means of the mean ct values in the housekeeping gene
ref_val <- Gmean_ct_target -Gmean_ct_hk     #calculate the reference ∆ct value
ref_val
hd1$delta_ct <- (hd1$target_mean_ct - hd1$hk_mean_ct) #create new column for  ∆ct to contain the newly calculated values 
hd1$delta_delta_ct <- (hd1$delta_ct - ref_val)  #create new column for  ∆∆ct to contain the newly calculated values 
hd1$expression <- (2^-hd1$delta_delta_ct) #calculate the expression data 
head(hd1)

#calculate the mean and standard errors (SE)
tab1_se<-summarise(
  group_by(hd1, tissue), 
  mean = mean(expression), 
  SE =(sd(expression))/sqrt(4)) #multiple statistics can be calculated within summarise 
tab1_se

#calculate relative expression: 
tab2_exp<-summarise(
  group_by(hd1, tissue), 
  target_mean = mean(target_mean_ct), 
  hk_mean = mean(hk_mean_ct)
)
tab2_exp

tab2_exp$delta_ct <- (tab2_exp$target_mean- tab2_exp$hk_mean) #add column to the tab2_exp dataframe for the ∆CT value
tab2_exp

ref_tab <- filter(tab2_exp, tissue == "root")
ref_val <- ref_tab$delta_ct
ref_val



tab2_exp$delta_delta_ct <- (tab2_exp$delta_ct - ref_val ) #add column to the tab2_exp dataframe for the ∆∆CT value
tab2_exp$expression <- (2^(-tab2_exp$delta_delta_ct)) #add column to the tab2_exp dataframe for the relative expression value
tab2_exp

complete_table <- data.frame(tab2_exp, tab1_se$SE)
complete_table


ggtitle.1 <- expression(paste(bolditalic("HDG5")))

hd1.plot <-ggplot(complete_table, aes(x = tissue, y = expression)) + 
  geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = tissue, ymin = expression - se, ymax = expression + se, width=0.1), size = 1) +
  xlab("") + ylab("Relative expression") + ggtitle(ggtitle.1) +
  scale_x_discrete(limits = c("root", "leaf", "flower"), 
                   labels = c("Root","Leaf", "Flower")) + 
  scale_y_continuous(
    limits = c(0,32),
    breaks=c(0,4,8,12,16,20,24,28,32)) + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  annotate("text", x =2, y = 5, label = "*", size = 8) +
  annotate("text", x = 3, y = 28, label = "***", size = 8) 

hd1.plot 

ggsave("hdg5_tissuetype.jpeg", device = "jpeg", units="in", width=5, height=4, dpi=300)

#############################################################################################

### This dataset includes 2 of 3 replicates for the trichome tissue subset ###
### Because there are only 2 reps, standard error and other stats cannot be calculated ###


urlfile1e="https://raw.githubusercontent.com/apicellap/data/main/gang_tissue_data_raw.csv"
gang_raw<-read.csv(url(urlfile1e)) 
gang_raw$Sample <- as.factor(gang_raw$Sample)

#qPCR data analysis: 
ref_df<-filter(gang_raw, tissue == "trichome") #subset tissue by root to isolate and do calculations on its mean ct values 
Gmean_ct_target <- mean(ref_df$target_mean_ct) #calculates means of the mean ct values in the target gene 
Gmean_ct_hk <- mean(ref_df$hk_mean_ct)   #calculates means of the mean ct values in the housekeeping gene 
ref_val <- Gmean_ct_target -Gmean_ct_hk     #calculate the reference ∆ct value
ref_val #print value in console 
gang_raw$delta_ct <- (gang_raw$target_mean_ct - gang_raw$hk_mean_ct) #create new column for  ∆ct to contain the newly calculated values 
gang_raw$delta_delta_ct <- (gang_raw$delta_ct - ref_val)  #create new column for  ∆∆ct to contain the newly calculated values 
gang_raw$expression <- (2^-gang_raw$delta_delta_ct) #calculate the expression data 


tab1_se_gn1<-summarise(
  group_by(gang_raw, tissue), 
  mean = mean(expression), 
  SE =(sd(expression))/sqrt(n())) # n() counts the number of observations in each group 
tab1_se_gn1


tab2_exp_gn1<-summarise(
  group_by(gang_raw, tissue), 
  target_mean = mean(target_mean_ct), 
  hk_mean = mean(hk_mean_ct))
tab2_exp_gn1

tab2_exp_gn1$delta_ct <- (tab2_exp_gn1$target_mean- tab2_exp_gn1$hk_mean)
tab2_exp_gn1

ref_tab_gn1 <- filter(tab2_exp_gn1, tissue == "trichome")
ref_val_gn1 <- ref_tab_gn1$delta_ct
ref_val_gn1


tab2_exp_gn1$delta_delta_ct <- (tab2_exp_gn1$delta_ct - ref_val_gn1 )
tab2_exp_gn1$expression <- (2^(-tab2_exp_gn1$delta_delta_ct))
tab2_exp_gn1

complete_tab_gn1 <- data.frame(tab2_exp_gn1, tab1_se_gn1$SE)
complete_tab_gn1

complete_tab_gn1 = rename(complete_tab_gn1, se = tab1_se_gn1.SE)
complete_tab_gn1

gang.plot <-ggplot(complete_tab_gn1, aes(x = tissue, y = expression)) + geom_bar(stat = 'identity', fill = "#82807E", color = "black", size =1) + 
  geom_errorbar(aes(x = tissue, ymin = expression - se, ymax = expression + se, width=0.1), size = 1) +
  xlab("") + ylab("") + ggtitle(ggtitle.1) +
  scale_x_discrete(limits = c("trichome", "trichomeless", "flower"), 
                   labels = c("Trichome","Trichomeless \n flower", "Flower")) + 
  scale_y_continuous(
    limits = c(0,600),
    breaks=seq(0,600, by = 60)) + 
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
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) 
gang.plot 

#############################################################################################

tis.plot<-(hd1.plot|gang.plot)
tis.plot

ggsave("02_figure.jpeg", device = "jpeg", units="in", width=12, height=6, dpi=300)


