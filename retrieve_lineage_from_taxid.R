#!/usr/bin/env Rscript
#retrieve lineage information for NCBI taxids created with diamond outfmt 102
#written by Philipp Resl, Nov 2017

#install.packages("taxize")
require(taxize)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("Need a file to work with.")
}

data <- read.csv(file = args[1], sep="\t",header=FALSE)
ids <- data$V2
for (i in 1:length(ids)) { #deal with multiple taxon ids, always select first!
	if (grepl(';',ids[i],fixed=TRUE)==TRUE) {
		ids[i]<-unlist(strsplit(as.character(ids[i]), split=";"))[1]
		print("Multiple taxonids found, will use first entry")
        }
}
data$V2 <- ids
ids <- as.character(unique(ids))
class <- classification(ids, db="ncbi")
data["V4"] <- NA 

for (i in 1:nrow(data)) {
    if (data[i,2] %in% ids) {
      if ("name" %in% colnames(class[[toString(data[i,2])]])){
        data[i,4] <-paste(class[[toString(data[i,2])]]$name, collapse=" ")
      }
      else {
        data[i,4] <- NA
      }
    }
}
write.table(file="",data,sep = "\t", quote=F, row.names=F, col.names=F)
#cat(data, row.names=F, col.names=F, sep="\t")
