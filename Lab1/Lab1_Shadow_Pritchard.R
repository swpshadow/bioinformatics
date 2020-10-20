#A
#1
seq(1,33,2)

#2
seq(7,40, length.out = 15)

#3
my.dna <- sample(c("A", "C", "G", "T"), size= 20, replace=TRUE)

#4
length(which(my.dna == "A"))

#5
barplot(table(my.dna), xlab = "DNA Bases", ylab = "count")
pie(table(my.dna))


#C
my.dna <- sample(c("A", "C", "G", "T"), size= 20, replace=TRUE, prob=c(.1,.4,.4,.1))
#only need to set working directory on my PC. 
#setwd("/home/file_path")
apoe.vec <- fasta2vec("Lab1/apoe.fasta")
#library(seqinr)


#function to read fasta file into a vector. 
fasta2vec <- function(fasta.file){
  my.fasta <- read.fasta(file=fasta.file, as.string=TRUE)
  my.fasta.string <- my.fasta[[1]][1]
  my.fasta.list <- strsplit(my.fasta.string,"")
  my.fasta.vec <- unlist(my.fasta.list)
}
#3
length(apoe.vec)
#4
apoe.vec[1:20]

#5 prints counts of each nucleotide
for (x in c('a', 'c', 'g', 't')) {
  sprintf("%s: %d", x, length(which(apoe.vec == x)))
}


#code to print the frequencies of the nucs in apoe.vec
for (x in c('a', 'c', 'g', 't')) { 
  print(paste0(x, ": ", length(which(apoe.vec == x))/length(apoe.vec)))
}

#6 makes bar and pie charts for apoe.vec
pie(table(apoe.vec))

barplot(table(apoe.vec), xlab = "DNA Bases", ylab = "Count")


#D
#1 GC count
(sum(apoe.vec=='g') + sum(apoe.vec=='c'))/length(apoe.vec)
GC(apoe.vec)
#3
pair.counts<-seqinr::count(apoe.vec,2)
#4
barplot(pair.counts)

#E loads the new file and prints the probabilities of each nuc
dna.vec <- fasta2vec("Lab1/DNA.fasta")
for (x in c('a', 'c', 'g', 't')) { 
    print(paste0(x, ": ", length(which(dna.vec == x))/length(dna.vec)))
}