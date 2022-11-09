#fungal proteome comparison from made proteome
#instal packages below before librarying them
library(readxl)
library(factoextra)
library(missMDA)
library(corrplot)
library(ggfortify)
library(BayesFactor)
library(ggpubr)
library(plyr)
library(stringr)


#import the genome sheet 
Fungal <- read_excel("Data/Fungal.EXCEL/try1build.xlsx",sheet = "Genomes")#, col_types = c("text","numeric", "numeric", "numeric", "numeric",
                                                                           #           "numeric", "numeric","numeric", "numeric", "numeric", 

#read in orthologs from the excel file

Orthologs <- read_excel("Data/Fungal.EXCEL/try1build.xlsx",sheet = "Orthologs")




#data engineering for the orthologs dataset
Trimmed<-as.data.frame(Orthologs[,-c(1:4)])
rownames(Trimmed)<-Orthologs$...1


#create presence absence matrix of reactions
Trimmed[!is.na(Trimmed)] <- 1
Trimmed[is.na(Trimmed)] <- 0

#make columns numeric
Trimmed<-apply(Trimmed, 2, as.numeric)   
#transpose matrix for analysis 
Trimmed<-t(Trimmed)






#create PCA data to use in the following visualizations
ORTHOPCA<-prcomp(Trimmed)

#assess percent variance explained by each eigenvalue (principal component)


fviz_eig(ORTHOPCA, choice = "variance", addlabels = T)


#visualize PCA of all fungi

PCA.IND <-
  fviz_pca_ind(
    ORTHOPCA,
    axes = c(1, 2),
    geom.ind = "point",
    pointshape = 19,
    pointsize = 3,
    #label = SAMMETA$CORE12,
   col.ind = rownames(Trimmed),
    mean.point = F
  )

PCA.IND


#ignome what is below here for now and proceed to ortholog_geneCounts

#compare Ortholog data 

#create a data set
Fungal<-as.data.frame(Fungal)

Fungal<-rowname.true(Fungal)
#for loop to create a matrix for a heat map
for (I in 1:ncol(Fungal)) {
  print(I)
  Fungal[,I]<-Fungal[,I]/Fungal[I,I]
  
}


