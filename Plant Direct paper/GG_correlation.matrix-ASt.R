#install/load packages:
install.packages("corrplot")
install.packages("Hmisc")
library(corrplot)
library(Hmisc)

#Download the data: 
urlfile4="https://raw.githubusercontent.com/apicellap/data/main/CtPharma_age_corr.csv"
GG_mat<-read.csv(url(urlfile4)) #data from CtPharma Gorilla Glue (gene expression & metabolite levels)

#Dataframe modifications: 
names(GG_mat)[names(GG_mat) == 'THC'] <- 'perc.THC' #changes name of variable
names(GG_mat)[names(GG_mat) == 'CBG'] <- 'perc.CBG'
GG_mat$dec.THC <- (GG_mat$perc.THC/100) #creates new column where percentage is converted to decimal, required by arcsin() calculation
GG_mat$dec.CBG <- (GG_mat$perc.CBG/100)
GG_mat$THC <-asin(sqrt(GG_mat$dec.THC)) #arcsin calculation. 
GG_mat$CBG<-asin(sqrt(GG_mat$dec.CBG))



#Create the correlation matrix: 
x <- subset(GG_mat, select = c(2:9, 14:15)) #subset the GG_mat dataframe into a new object called 'x'
GGcor <-rcorr(as.matrix(x)) #coerce x into being a matrix
GGcor #print the matrix
GGcor$P[is.na(GGcor$P)]<-1 #replace missing values with 1s 

#this is the correct order of making the plot
#you have to open the file, then run the plot, and then close the file with dev.off
tiff(file="GG_arcsin-matrix.tiff",
     width=6, height=4, units="in", res=200)

corrplot(GGcor$r, #strength of correlation (info for size and color of the shape)
         method="circle", #each correlation represented by a shape, in this case a circle 
         type="upper", #orientation of plot (upper portion of the diagnonal)
         p.mat=GGcor$P, #p values for correlation, is the correlation a statistically significant, if so, then how significant? 
         insig="label_sig", 
         pch.col="white",  #color of symbol 
         sig.level = c(.001,.01,.05), #p value cutoff points 
         pch.cex=1.25, #size of symbol
         tl.srt = 45) #tilt labels above graph at a 45˚ angle 
dev.off()

ggsave("corr-matrix.tiff", units="in", width=5, height=5, dpi=175, compression = 'lzw')

