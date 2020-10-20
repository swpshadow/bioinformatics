# load gene expression data


#setwd("/home/swpshadow/Documents/Bio/Lab6")
load("sense.filtered.cpm.Rdata")  # setwd!

# load phenotype (mdd/hc) data
subject.attrs <- read.csv("Demographic_symptom.csv", 
                          stringsAsFactors = FALSE)
library(dplyr)
# grab intersecting X (subject ids) and Diag (Diagnosis) from columns
phenos.df <- subject.attrs %>% 
  filter(X %in% colnames(sense.filtered.cpm)) %>%
  dplyr::select(X, Diag)  
mddPheno <- as.factor(phenos.df$Diag)  

# Normalized and transform
library(preprocessCore)
mddExprData_quantile <- normalize.quantiles(sense.filtered.cpm)
mddExprData_quantileLog2 <- log2(mddExprData_quantile)
# attach phenotype names and gene names to data 
colnames(mddExprData_quantileLog2) <- mddPheno  
rownames(mddExprData_quantileLog2) <- rownames(sense.filtered.cpm)

# coefficient of variation filter sd(x)/abs(mean(x))
CoV_values <- apply(mddExprData_quantileLog2,1,
                    function(x) {sd(x)/abs(mean(x))})
# smaller threshold, the higher the experimental effect relative to the 
# measurement precision
thresh <- .02
sum(CoV_values<.thresh)
# there is one gene that has 0 variation -- remove
sd_values <- apply(mddExprData_quantileLog2,1, function(x) {sd(x)})
rownames(mddExprData_quantileLog2)[sd_values==0]  
# filter the data matrix 
GxS.covfilter <- mddExprData_quantileLog2[CoV_values< thresh & sd_values>0,]
dim(GxS.covfilter)

# convert phenotype to factor
pheno.factor <- as.factor(colnames(GxS.covfilter))
pheno.factor
str(pheno.factor)
levels(pheno.factor)


############## A ##################
# 1.ex.data$map
snprow <- which(ex.data$map$snp.names=="rs630969")
ex.data$map$allele.1[snprow]
ex.data$map$allele.2[snprow]

mddCorr<-cor(t(GxS.covfilter))  # correlation between genes
thresh<-.7 # controls sparsity of network
# threshold and turn T/F to 1/0
adjMat <- (abs(mddCorr)>thresh)+0
diag(adjMat) <- 0  # remove self-connections
rownames(adjMat) <- row.names(GxS.covfilter)
colnames(adjMat) <- row.names(GxS.covfilter)

# 2.
library(igraph)
ig <- graph.adjacency(adjMat, mode = "undirected")
plot(ig, vertex.size=1, vertex.label.color = "black", 
     edge.width=2, vertex.label=NA)
igDegrees <- rowSums(adjMat) # degree vector for nodes
hist(igDegrees)

# 3.
# remove genes that have no connections
adjMat.connected <- adjMat[igDegrees!=0,igDegrees!=0]
ig.connected <- graph.adjacency(adjMat.connected, mode = "undirected")
plot(ig.connected, vertex.size=1, vertex.label.color = "black", 
     edge.width=2, vertex.label=NA)
hist(rowSums(adjMat.connected))

############## B ##################
# 4.

# fast community detection algorithm
igc.clusts <- fastgreedy.community(ig.connected)

length(igc.clusts)
sizes(igc.clusts)
igc.membership <- membership(igc.clusts)
length(igc.membership)


#5.

igc.colors <- igc.membership
color.pallet <- rainbow(length(igc.clusts)) # discrete color for each cluster
for (i in (1:length(igc.clusts))){  # change membership to the color pallet
  igc.colors[i==igc.colors]<-color.pallet[i]
}
plot(ig.connected, vertex.size=1, vertex.label.color = igc.membership,
     vertex.color = igc.membership, edge.width=2)

# 6.
clust1.genes<-names(igc.membership)[igc.membership==1]

write.table(clust1.genes,file="clust1.txt",row.names=F,col.names=F,quote=F) 

msigdb = read.table("clust1.txt")

############## C ##################

# 7. 8.
library(snpStats) # install first

ex.data <- read.pedfile(file="extra.ped", snps="extra.map")
ex.data$fam
phenotype <- ex.data$fam$affected-1  # change pheno from 1/2 to 0/1
genotypes <- ex.data$genotypes       # encoded as AA/AB/BB
print(sum(phenotype==1))


# 9. 
snp.ids <- as.character(ex.data$map$snp.names)

genotypes.df <- data.frame(as(genotypes, "character"))
table(genotypes.df$rs630969)
table(phenotype,genotypes.df$rs630969)

# 10.
ex.data$map
snprow <- which(ex.data$map$snp.names=="rs630969")
ex.data$map$allele.1[snprow]
ex.data$map$allele.2[snprow]
