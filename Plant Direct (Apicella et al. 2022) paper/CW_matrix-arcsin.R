#install/load packages:
install.packages("corrplot")
install.packages("Hmisc")
library(corrplot)
library(Hmisc)

urlfile2a="https://raw.githubusercontent.com/apicellap/data/main/CherryWine_correlation.csv"
CW_mat<-read.csv(url(urlfile2a))

#Dataframe modifications: 
names(CW_mat)[names(CW_mat) == 'THC'] <- 'perc.THC' #changes name of variable
names(CW_mat)[names(CW_mat) == 'CBC'] <- 'perc.CBC' 
names(CW_mat)[names(CW_mat) == 'CBG'] <- 'perc.CBG' 
names(CW_mat)[names(CW_mat) == 'CBD'] <- 'perc.CBD'
CW_mat$dec.THC <- (CW_mat$perc.THC/100) #creates new column where percentage is converted to decimal, required by asin() calculation
CW_mat$dec.CBC <- (CW_mat$perc.CBC/100)
CW_mat$dec.CBG <- (CW_mat$perc.CBG/100)
CW_mat$dec.CBD <- (CW_mat$perc.CBD/100)
CW_mat$THC <-asin(sqrt(CW_mat$dec.THC)) #arcsin calculation. 
CW_mat$CBC <-asin(sqrt(CW_mat$dec.CBC))
CW_mat$CBG <-asin(sqrt(CW_mat$dec.CBG))
CW_mat$CBD <-asin(sqrt(CW_mat$dec.CBD))

#Create the correlation matrix: 
y <- subset(CW_mat, select = c(2:9, 19:22)) #subset the GG_mat dataframe into a new object called 'y'
cor4 <-rcorr(as.matrix(y)) #coerce y into being a matrix
cor4 #print the matrix
cor4$P[is.na(cor4$P)]<-1 #replace missing values with 1s 

y1 <- subset(CW_mat) #subset the GG_mat dataframe into a new object called 'y'
cor5 <-rcorr(as.matrix(y1)) #coerce y into being a matrix
cor5 #print the matrix
cor5$P[is.na(cor5$P)]<-1 #replace missing values with 1s 

#this is the correct order of saving the plot to a directory 
  #you have to open the file using the tiff() function, 
  #then run the plot, and then close the file with dev.off
tiff(file="CW_matrix-arcsin.tiff",
     width=6, height=4, units="in", 
     res=310, compression = "lzw") 

corrplot(cor4$r, #strength of correlation (info for size and color of the shape)
         method="circle", #each correlation represented by a shape, in this case a circle 
         type="upper", #orientation of plot (upper portion of the diagnonal)
         p.mat=cor4$P, #p values for correlation, is the correlation a statistically significant, if so, then how significant? 
         insig="label_sig", 
         pch.col="white",  #color of symbol 
         sig.level = c(.001,.01,.05), #p value cutoff points 
         pch.cex=1.25, #size of symbol
         tl.srt = 45) #tilt labels above graph at a 45˚ angle 
dev.off()
