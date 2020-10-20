#############   A   #####################
# load gene expression data
setwd("/home/swpshadow/Documents/Bio/Lab4")
load("sense.filtered.cpm.Rdata")

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

dim(mddExprData_quantileLog2)

#############   B   #####################
# coefficient of variation filter sd(x)/abs(mean(x))
CoV_values <- apply(mddExprData_quantileLog2,1,
                    function(x) {sd(x)/abs(mean(x))})
# smaller threshold, the higher the experimental effect relative to the 
# measurement precision
sum(CoV_values<.045)
# there is one gene that has 0 variation -- remove
sd_values <- apply(mddExprData_quantileLog2,1, function(x) {sd(x)})
rownames(mddExprData_quantileLog2)[sd_values==0]  
# filter the data matrix 
GxS.covfilter <- mddExprData_quantileLog2[CoV_values<.045 & sd_values>0,]
dim(GxS.covfilter)


#############   C   #####################
# convert phenotype
pheno.factor <- as.factor(colnames(GxS.covfilter))
pheno.factor
str(pheno.factor)
levels(pheno.factor)

myrow <- 2  # first pick a gene row index to test
mygene<-rownames(GxS.covfilter)[myrow]
mygene

# a. traditional R interface 
mdd <- GxS.covfilter[myrow,pheno.factor=="MDD"]
hc <- GxS.covfilter[myrow,pheno.factor=="HC"]
t.result <- t.test(mdd,hc)
t.result  

# b. formula interface ~ saves a step
t.result <- t.test(GxS.covfilter[myrow,] ~ pheno.factor)
t.result
t.result$p.value
t.result$statistic

#5. plot the data
library(ggplot2)
# create data frame for gene
mygene.data.df <- data.frame(gene=GxS.covfilter[myrow,],phenotype=pheno.factor)
# boxplot
p <- ggplot(mygene.data.df, aes(x=phenotype, y=gene, fill=phenotype)) + stat_boxplot(geom ='errorbar') + geom_boxplot()
p <- p + xlab("MDD versus HC") + ylab(mygene)
p


#############   D   #####################

# Put it all together into a function to run in loop. 
# First write a function that computes t-test for one gene.
# i is the data row for the gene
ttest_fn <- function(i){
  mygene <- rownames(GxS.covfilter)[i]
  t.result <- t.test(GxS.covfilter[i,] ~ pheno.factor)
  tstat <- t.result$statistic
  pval <- t.result$p.value
  # return vector of three things for each gene    
  c(mygene, tstat, pval) 
} 

ttest_fn(2) 


# initialize an empty matrix to store the results
ttest_allgene.mat <- matrix(0,nrow=nrow(GxS.covfilter), ncol=3)
# run analysis on all gene rows
for (i in 1:nrow(GxS.covfilter)){
  ttest_allgene.mat[i,] <- ttest_fn(i) 
}
# convert matrix to data frame and colnames
ttest_allgene.df <- data.frame(ttest_allgene.mat)
colnames(ttest_allgene.df) <- c("gene ", "t.stat", "p.val")

library(dplyr)          # sort based on p-value
ttest_allgene.sorted <- ttest_allgene.df                 %>% 
  mutate_at("p.val", as.character) %>%
  mutate_at("p.val", as.numeric)   %>%
  arrange(p.val) # sort 
ttest_allgene.sorted[1:10,] # look at top 10

#8. find data row index of top gene name
myrow <- which(ttest_allgene.df$gene=="MDGA1")
mygene<-rownames(GxS.covfilter)[myrow]

mygene.data.df <- data.frame(gene=GxS.covfilter[myrow,],phenotype=pheno.factor)
# boxplot
p <- ggplot(mygene.data.df, aes(x=phenotype, y=gene, fill=phenotype)) + stat_boxplot(geom ='errorbar') + geom_boxplot()
p <- p + xlab("MDD versus HC") + ylab(mygene)
p


#9

top_cutoff <- 200
top_genes <- as.character(ttest_allgene.sorted[1:top_cutoff,1])
write.table(top_genes, sep="\t", file="", quote=F, row.names=F, col.names=F)



