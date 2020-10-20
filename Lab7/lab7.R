library(snpStats) #BiocManager::install("snpStats") 

#setwd("/home/swpshadow/Documents/Bio/Lab7")

ex.data <- read.pedfile(file="extra.ped", snps="extra.map")
ex.data$fam
phenotype <- ex.data$fam$affected-1  # change pheno from 1/2 to 0/1
genotypes <- ex.data$genotypes       # encoded as AA/AB/BB
snp.ids <- as.character(ex.data$map$snp.names)
colnames(genotypes.df) <- snp.ids
genotypes.df <- data.frame(as(genotypes, "character"))
# observed contingency table for SNP rs630969
table(phenotype,genotypes.df$rs630969,
      dnn=c("phenotype","genotype")) # dnn dimension names of table 


### A ###

# creates list of observed contingency tables for all SNPs
# sapply acts on each column of genotypes.df
observed.tables.list <- sapply(genotypes.df, function(x) table(phenotype,x,dnn=c("phenotype","genotype")))


test.table <- observed.tables.list$rs634228 # grab one table
genoMarg.vec  <-  colSums(test.table)                     # margin vector
phenoMarg.vec <- rowSums(test.table)                     # margin vector        
totalSubj <- sum(genoMarg.vec)                     # total subjects

expect.test <- outer(phenoMarg.vec,genoMarg.vec/totalSubj,'*')


# Fisher exact test (chi-square test) for all SNPs
fish_fn <- function(i){  
  cbind(snp.ids[i], fisher.test(observed.tables.list[[i]])$p.value)    
}

# apply fisher exact test to all SNPs
fish.df <- data.frame(t(sapply(1:ncol(genotypes.df), fish_fn)))
colnames(fish.df) <- c("rs", "p_value")

# sort SNPs by Fisher exact p-value
library(dplyr)
fish.results <- fish.df                            %>% 
  mutate_at("p_value", as.character) %>%
  mutate_at("p_value", as.numeric)   %>%
  arrange(p_value) 
print(fish.results)



### B ###
library(ggplot2)
i<-8
A1<-ex.data$map$allele.1[i]
A2<-ex.data$map$allele.2[i]
geno.labels <- c(paste("A","A",sep="/"),paste("A","B",sep="/"),paste("A","B",sep="/"))

# data from the one SNP
oneSNP.df <- data.frame(cbind(genotypes.df[[i]],as.numeric(phenotype)))
colnames(oneSNP.df) <- c("genotypes","probability")

lr.plot <- ggplot(oneSNP.df, aes(x=genotypes, y=probability)) + 
  geom_point(position = position_jitter(w = 0.1, h = 0.1)) +
  stat_smooth(method="glm", method.args = list(family = "binomial")) +
  xlim(geno.labels) + ggtitle(snp.ids[i])
print(lr.plot)


library(broom) # for tidy function. make sure installed
pheno.factor <- factor(phenotype,labels=c(0,1))
i<-8
lr <- glm(pheno.factor~genotypes.df[[i]],family=binomial)
td.lr <- tidy(lr) 

pval_vec <- td.lr$p.value # vector of $p.values from td.lr

coef_vec <- td.lr$estimate # vector of $estimates

cbind(snp.ids[i], coef_vec[1], coef_vec[2], coef_vec[3], pval_vec[1], pval_vec[2], pval_vec[3])  



LR.fn <- function(i){
  
  lr <- glm(pheno.factor~genotypes.df[[i]],family=binomial)
  
  td.lr <- tidy(lr)
  
  pval_vec <- td.lr$p.value # vector of $p.values from td.lr
  
  coef_vec <- td.lr$estimate # vector of $estimates
  
  cbind(snp.ids[i], coef_vec[1], coef_vec[2], coef_vec[3], pval_vec[1], pval_vec[2], pval_vec[3])  
} 

# apply Logistic Regression model to all SNPs
LRresults.df <- data.frame(t(sapply(1:ncol(genotypes.df), LR.fn)))

# add column names to results data frame
colnames(LRresults.df) <- c("rs", "AAintercept", "ABcoef", "BBcoef", "AA.pval", "AB.pval", "BB.pval")

# sort LR results by the p-value of the BB homozygous coefficient
# tidy made $p_value a factor and when you try to convert directly to numeric
# as.numeric turns factors into integer and this messes up sorting
# especially scientific notation

lr.results.sorted <- LRresults.df %>% 
  mutate_at("BB.pval", as.character) %>%
  mutate_at("BB.pval", as.numeric) %>%
  arrange(BB.pval) 

as.matrix(lr.results.sorted %>% pull(rs,BB.pval))
