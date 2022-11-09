#using orthologs to look and gene count
library(qgraph)
library(vegan)
#create copy of orthologs to store results

Orthologs2<-Orthologs
#for loop to convert orthologs to frequency data 


#subsets each row
for(i in 1:nrow(Orthologs)){
  #in each row go through each column 
  for(j in 5:ncol(Orthologs)) {
    #for each cell we are replace genes with a number, 
    #for ease we are replacing them with a character number
    #each gene is speparated by a ; so count them
    #if it is NA  replace with "0"
    if(is.na(str_count(Orthologs[2,18], ";")))
      Orthologs2[i,j]<-"0"
    #if there are ; add number of ; plus one for the number of orthologs
    else{
    Orthologs2[i,j]<-as.character(str_count(Orthologs[i,j], ";")+1)
    }
    
  }
  
  
}

#convert character columns to numeric 

Orthologs2[,5:ncol(Orthologs2)]<-as.data.frame(lapply(Orthologs2[,5:ncol(Orthologs2)],as.numeric))
  
#replace NA with 0 

Orthologs2[,5:ncol(Orthologs2)][is.na(Orthologs2[,5:ncol(Orthologs2)])] <- 0

#change tibble to data frame
Orthologs2<-as.data.frame(Orthologs2)

#remove extra data and change row names to functions
row.names(Orthologs2)<-make.names(Orthologs2$`representative function`, unique = TRUE)

Orthologs2<-Orthologs2[,-c(1:4)]

View(head(Orthologs2))

#transpose so genome is the row 

T_Orthologs2<-t(Orthologs2)

#create PCA data to use in the following visualizations
ORTHOPCA<-prcomp(T_Orthologs2)

#assess precent variance explained by each eigenvalue (principal component)


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


#make a distance matrix Eucledian for now

Orthologs2_dist_Euc <- dist(T_Orthologs2, method="euclidean")


plot(
  hclust(Orthologs2_dist_Euc),
  hang = -1,
  main = "Species clustered by Jaccard similarity",
  axes = FALSE, ylab = ""
)



#diversity measures
library(vegan)


T_PAV<-(T_Orthologs2)

#inverse simpson is a measure of alpha diversity
PAV.invsimpson<-diversity(T_PAV, index = "invsimpson")
#

#function that makes a plot
make.sorted.plot <- function(x){
  ordered <- sort(x, T)
  plot(
    ordered,
    col = terrain.colors(10),
    xaxt = "n", pch = 16, cex = 2,
    ylim = c(min(ordered)*0.5, max(ordered)),
    xlim = c(0, length(x)+1),
    ylab = "Diversity measure", xlab = "Sites",
    main = substitute(x))
  text(ordered,
       names(ordered),
       srt = -75,
       pos = 4) }

#plot; the higher the measure the greater the number of effective groups. 
#e.g. a low value = among orthologs functional groups this species has a low number of groups relative to whats in the full data set
#for instance 200 total ortholog groups, if the value is low the organism might have orthologs pertaining to a small number of the groups
#the opposite is true where a high number means they have orthologs in every group
#this also has an assumption of evenness among species 
#which is true in this instance. e.g. if we take any of the other types of diversity indicies the values are basically the sam e
#its interesting to B. porri is so low 
make.sorted.plot(PAV.invsimpson)


PAV.jaccard <- vegdist(T_PAV, method = "jaccard",binary = F)

#here is a dendrograme using jaccard similarity
#calculated by dividing the number of observations of orthologs in two species 
#by the number of observations of orthologs in either species
#species closer together have a higher number of orthologs in the functional groups 
#when standardized to the number of orthologs in each species
#branch length corresponds to differences of the node
#e.g. short branch = very similar, long branch = very different 
# Dissimilarity dendrogram
dune.dist<-vegdist(T_Orthologs2, method = "jaccard")
dune.hc<-hclust(dune.dist)
plot(dune.hc)



#lets make a heat map of similarities


#lets make a heat map of dissimilarities
heatmap(as.matrix(dune.dist),symm = T)

#lets make a network map map of similarities
plot(dune.dist)
dist_m <- as.matrix(dune.dist)
dist_mi <- 1-dist_m # minus 1, as qgraph takes similarity matrices as input

jpeg('Botrytis_Orthologs_forcedraw.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring',tuning=0.5, vsize=5, node.width=2, labels=T)
dev.off()
#the numbers on the graph correspond to the order of the rows or columns(its symmetric) in 
View(dist_m)



##### we are gonna do the same analyses minus the outgroups

View(head(Orthologs2))

Orthologs3<-Orthologs2[,-c(11,14:16)]

Orthologs3<-Orthologs3[rowSums(Orthologs3[])>0,]


#remove all rows = to zero 

T_Orthologs3<-t(Orthologs3)

#create PCA data to use in the following visualizations
ORTHOPCA<-prcomp(T_Orthologs3)

#assess precent variance explained by each eigenvalue (principal component)


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
    col.ind = rownames(T_Orthologs3),
    mean.point = F
  )

PCA.IND


#make a distance matrix Eucledian for now

Orthologs3_dist_Euc <- dist(T_Orthologs3, method="euclidean")


plot(
  hclust(Orthologs3_dist_Euc),
  hang = -1,
  main = "Species clustered by Jaccard similarity",
  axes = FALSE, ylab = ""
)



#diversity measures
library(vegan)


T_PAV<-(T_Orthologs3)

#inverse simpson is a measure of alpha diversity
PAV.invsimpson<-diversity(T_PAV, index = "invsimpson")
#

#function that makes a plot
make.sorted.plot <- function(x){
  ordered <- sort(x, T)
  plot(
    ordered,
    col = terrain.colors(10),
    xaxt = "n", pch = 16, cex = 2,
    ylim = c(min(ordered)*0.5, max(ordered)),
    xlim = c(0, length(x)+1),
    ylab = "Diversity measure", xlab = "Sites",
    main = substitute(x))
  text(ordered,
       names(ordered),
       srt = -75,
       pos = 4) }

#plot; the higher the measure the greater the number of effective groups. 
#e.g. a low value = among orthologs functional groups this species has a low number of groups relative to whats in the full data set
#for instance 200 total ortholog groups, if the value is low the organism might have orthologs pertaining to a small number of the groups
#the opposite is true where a high number means they have orthologs in every group
#this also has an assumption of evenness among species 
#which is true in this instance. e.g. if we take any of the other types of diversity indicies the values are basically the sam e
#its interesting to B. porri is so low 
make.sorted.plot(PAV.invsimpson)


PAV.jaccard <- vegdist(T_PAV, method = "jaccard",binary = F)

#here is a dendrograme using jaccard similarity
#calculated by dividing the number of observations of orthologs in two species 
#by the number of observations of orthologs in either species
#species closer together have a higher number of orthologs in the functional groups 
#when standardized to the number of orthologs in each species
#branch length corresponds to differences of the node
#e.g. short branch = very similar, long branch = very different 
# Dissimilarity dendrogram
dune.dist<-vegdist(T_Orthologs3, method = "jaccard")
dune.hc<-hclust(dune.dist)
plot(dune.hc)



#lets make a heat map of similarities


#lets make a heat map of dissimilarities
heatmap(as.matrix(dune.dist),symm = T)

#lets make a network map map of similarities
dist_m <- as.matrix(dune.dist)
 # minus 1, as qgraph takes similarity matrices as input
dist_mi <- 1-dist_m
jpeg('Botrytis_Orthologs_forcedraw.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring', vsize=5, node.width=2, labels=T)
dev.off()
#the numbers on the graph correspond to the order of the rows or columns(its symmetric) in 
View(dist_mi)
#in this graph fragrae and cinerea are the most differen tfrom one another 

##########adding one more to look at without porri



View(head(Orthologs2))

Orthologs3<-Orthologs2[,-c(14:16)]

Orthologs3<-Orthologs3[rowSums(Orthologs3[])>0,]


#remove all rows = to zero 

T_Orthologs3<-t(Orthologs3)

#create PCA data to use in the following visualizations
ORTHOPCA<-prcomp(T_Orthologs3)

#assess precent variance explained by each eigenvalue (principal component)


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
    col.ind = rownames(T_Orthologs3),
    mean.point = F
  )

PCA.IND


#make a distance matrix Eucledian for now

Orthologs3_dist_Euc <- dist(T_Orthologs3, method="euclidean")


plot(
  hclust(Orthologs3_dist_Euc),
  hang = -1,
  main = "Species clustered by Jaccard similarity",
  axes = FALSE, ylab = ""
)



#diversity measures



T_PAV<-(T_Orthologs3)

#inverse simpson is a measure of alpha diversity
PAV.invsimpson<-diversity(T_PAV, index = "invsimpson")
#

#function that makes a plot
make.sorted.plot <- function(x){
  ordered <- sort(x, T)
  plot(
    ordered,
    col = terrain.colors(10),
    xaxt = "n", pch = 16, cex = 2,
    ylim = c(min(ordered)*0.5, max(ordered)),
    xlim = c(0, length(x)+1),
    ylab = "Diversity measure", xlab = "Sites",
    main = substitute(x))
  text(ordered,
       names(ordered),
       srt = -75,
       pos = 4) }

#plot; the higher the measure the greater the number of effective groups. 
#e.g. a low value = among orthologs functional groups this species has a low number of groups relative to whats in the full data set
#for instance 200 total ortholog groups, if the value is low the organism might have orthologs pertaining to a small number of the groups
#the opposite is true where a high number means they have orthologs in every group
#this also has an assumption of evenness among species 
#which is true in this instance. e.g. if we take any of the other types of diversity indicies the values are basically the sam e
#its interesting to B. porri is so low 
make.sorted.plot(PAV.invsimpson)


PAV.jaccard <- vegdist(T_PAV, method = "jaccard",binary = F)

#here is a dendrograme using jaccard similarity
#calculated by dividing the number of observations of orthologs in two species 
#by the number of observations of orthologs in either species
#species closer together have a higher number of orthologs in the functional groups 
#when standardized to the number of orthologs in each species
#branch length corresponds to differences of the node
#e.g. short branch = very similar, long branch = very different 
# Dissimilarity dendrogram
dune.dist<-vegdist(T_Orthologs3, method = "jaccard")
dune.hc<-hclust(dune.dist)
plot(dune.hc)



#lets make a heat map of similarities


#lets make a heat map of dissimilarities
heatmap(as.matrix(dune.dist),symm = T)

#lets make a network map map of similarities
dist_m <- as.matrix(dune.dist)
# minus 1, as qgraph takes similarity matrices as input
dist_mi <- 1-dist_m
jpeg('Botrytis_Orthologs_forcedraw.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring', vsize=5, node.width=2, labels=T)
dev.off()
#the numbers on the graph correspond to the order of the rows or columns(its symmetric) in 
View(dist_mi)
#in this one aclada and hyacinthi are the most different



