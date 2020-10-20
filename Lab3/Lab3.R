setwd("/home/swpshadow/Documents/Bio")

load("./Lab3/sense.filtered.cpm.Rdata")
dim(sense.filtered.cpm) 
colnames(sense.filtered.cpm)


subject.attrs <- read.csv("./Lab3/Demographic_symptom.csv", stringsAsFactors = FALSE)
dim(subject.attrs)  # 160 subjects x 40 attributes
colnames(subject.attrs)  # interested in X (sample ids) and Diag (diagnosis)
subject.attrs$X
sum(subject.attrs$Diag == "MDD")

library(dplyr) # install.packages("dplyr")
# create a phenotype vector
# grab X (subject ids) and Diag (Diagnosis) from subject.attrs that 
# intersect %in% with the RNA-Seq data
phenos.df <- subject.attrs %>% 
  filter(X %in% colnames(sense.filtered.cpm)) %>%
  dplyr::select(X, Diag)  
colnames(phenos.df) # $Diag is mdd diagnosis
# grab Diag column and convert character to factor
mddPheno <- as.factor(phenos.df$Diag)  # this is our phenotype/class vector 

summary(mddPheno) # MDD -- major depressive disorder, HC -- healthy control


### raw cpm boxplots and histogram of one sample
boxplot(sense.filtered.cpm,range=0,ylab="raw probe intensity", main="Raw", names=mddPheno)

hist(sense.filtered.cpm[,1], freq=F, ylab="density", xlab="raw probe intensity", main="Raw Data Density for Sample 1")

### log2 transformed boxplots and histogram
boxplot(log2(sense.filtered.cpm), range=0,ylab="log2 intensity", main="Log2 Transformed", names=mddPheno)

hist(log2(sense.filtered.cpm[,1]), freq=F, ylab="density", xlab="log2 probe intensity", main="log2 Data Density for Sample 1")


# install quantile normalize
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("preprocessCore")
library(preprocessCore)

# apply quantile normalization
mddExprData_quantile <- normalize.quantiles(sense.filtered.cpm)

boxplot(mddExprData_quantile,range=0,ylab="raw intensity", main="Quantile Normalized")


mddExprData_quantileLog2 <- log2(mddExprData_quantile)
# add phenotype names to matrix
colnames(mddExprData_quantileLog2) <- mddPheno  

boxplot(mddExprData_quantileLog2,range=0,ylab="log2 intensity", main="Quantile Normalized Log2")

hist(log2(mddExprData_quantileLog2[,1]), freq=F, ylab="density", xlab="log2 probe intensity", main="log2 Quantile Normalized for Sample 1")

mean(mddExprData_quantileLog2[,1])
colMeans(mddExprData_quantileLog2)


# transpose data matrix and convert to data frame
# ggplot wants data frame and subjects as rows
expr_SxG <- data.frame(t(mddExprData_quantileLog2)) # Subject x Gene
colnames(expr_SxG) <- rownames(sense.filtered.cpm)  # add gene names

## MDS of subjects 
#d<-dist(expr_SxG)         # Euclidean metric
mddCorr<-cor(t(expr_SxG))  # distance based on correlation
d <- 1-mddCorr
mdd.mds <- cmdscale(d, k = 2)
x <- mdd.mds[,1]
y <- mdd.mds[,2]
mdd.mds.df <- data.frame(mdd.mds)
colnames(mdd.mds.df) <- c("dim1","dim2")
#install.packages("ggplot2") # if not already installed
#BiocManager::install("ggplot2")("ggplot2")
library(ggplot2)

p <- ggplot() # initialize empty ggplot object
p <- p + geom_point(data=mdd.mds.df, aes(x=dim1, y=dim2, color=mddPheno, shape=mddPheno), size=3)
p <- p + ggtitle("MDS") + xlab("Dim 1") + ylab("Dim 2")
print(p)


## hierarchical cluster of subjects
mddTree = hclust(as.dist(d))
mddTree$labels <- mddPheno
plot(mddTree)

num.clust <- 5
mddCuts <- cutree(mddTree,k=num.clust)
table(names(mddCuts),mddCuts)

#install.packages("dendextend")
library(dendextend)
dend <- as.dendrogram(mddTree)
dend=color_labels(dend, k=num.clust)
#dend <- dend %>% color_branches(k = 4) 
plot(dend) # selective coloring of branches AND labels 

library(BiocManager)
#BiocManager::install("WGCNA")
library(WGCNA)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
dynamicMods = cutreeDynamic(dendro = mddTree, distM = d,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 2)
mddColors = labels2colors(dynamicMods)
table(mddColors)
table(mddColors,names(mddCuts))
plotDendroAndColors(mddTree, mddColors, "Dynamic Clusters",
                    dendroLabels = NULL, # hang = -1,
                    addGuide = TRUE, #guideHang = 0.05,
                    main = "Clustering with WGCNA")
