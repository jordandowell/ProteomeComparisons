library(data.table)
library(reshape2)

#create distance figure for proteome comparison




#read in proteome comparison

protcomp_list <- read.delim("Data/protcomp_list.tsv", header=T)





Fungalgenes<-subset(protcomp_list, grepl('^[[:alpha:] ]+$', protcomp_list$Genome.name2))


indices <- which(protcomp_list$Genome.name2=='^[[:alpha:] ]+$')



X<-protcomp_list$Genome.name2
Y<-protcomp_list$Genome.name1
small<-protcomp_list[grep('^[[:alpha:] ]', X),, drop = FALSE]
smaller<-small[grep('^[[:alpha:] ]', Y),, drop = FALSE]
smallest<-smaller[!is.na(smaller$bit.score),]





#read in species names

Species<-read.csv("Data/4lettergenes.csv", header=T)
View(Species)



#add 4 letter column to smallest

smallest$fourletter<-substr(smallest$Genome.name2, start = 1, stop = 4)
View(head(smallest))


Proteome_species<-merge(smallest,Species, by.x = "fourletter", by.y = "Letter")


Proteome_species$Species<-as.factor(Proteome_species$Species)



#presence absenece variation


PAV<-table(Botrytis = Proteome_species$Genome.name1,Proteome_species$Species)

PAV<-as.data.frame(PAV)
#need to melt wide

PAV<-(dcast(PAV, Botrytis ~ Var2, value.var = "Freq"))

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
rowname.true <- function(df) {
  rownames(df) <- as.character(unlist(df[,1]))
  df[,-1]
}

PAV<-rowname.true(PAV)



#diversity measures
library(vegan)


T_PAV<-(t(PAV))


PAV.genes<-specnumber(T_PAV)

PAV.shannon<-diversity(T_PAV)
PAV.shannon.equitibility<-(PAV.shannon)/log(ncol(T_PAV))

View(PAV.shannon)
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


make.sorted.plot(PAV.shannon)
make.sorted.plot(PAV.genes)
make.sorted.plot(PAV.shannon.equitibility)


PAV.jaccard <- vegdist(T_PAV, method = "jaccard",binary = F)


plot(
  hclust(PAV.jaccard),
  hang = -1,
  main = "Species clustered by Jaccard similarity",
  axes = FALSE, ylab = ""
)

View(as.matrix(PAV.jaccard))
# Dissimilarity dendrogram
dune.dist<-vegdist(T_PAV)
dune.hc<-hclust(dune.dist)
plot(dune.hc)






####bbh ratio as a comparaible to percent cover 


View(Proteome_species)


PAV<-table(Botrytis = Proteome_species$Genome.name1,Proteome_species$Species)



BBH<-aggregate(Proteome_species$bbh.percent, by = list(Proteome_species$Genome.name1,Proteome_species$Species), FUN = max)





BBH<-(dcast(BBH, Group.1 ~ Group.2, value.var = "x"))

BBH[is.na(BBH)] <- 0



header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
rowname.true <- function(df) {
  rownames(df) <- as.character(unlist(df[,1]))
  df[,-1]
}

BBH<-rowname.true(BBH)


#remove everything but row maxima
library(matrixStats)
myvector <- rowMaxs(as.matrix(BBH))
BBH[BBH < myvector] <- 0

#diversity measures
library(vegan)


T_BBH<-(t(BBH))


BBH.genes<-specnumber(T_BBH)

BBH.shannon<-diversity(T_BBH)
BBH.shannon.equitibility<-(BBH.shannon)/log(ncol(T_BBH))


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


make.sorted.plot(BBH.shannon)
make.sorted.plot(BBH.genes)
make.sorted.plot(BBH.shannon.equitibility)


BBH.jaccard <- vegdist(T_BBH, method = "jaccard",binary = F)

heatmap(as.matrix(BBH.jaccard))

plot(
  hclust(BBH.jaccard),
  hang = -1,
  main = "Species clustered by Jaccard similarity",
  axes = FALSE, ylab = ""
)


# Dissimilarity dendrogram
dune.dist<-vegdist(T_BBH)
dune.hc<-hclust(dune.dist)
plot(dune.hc)





