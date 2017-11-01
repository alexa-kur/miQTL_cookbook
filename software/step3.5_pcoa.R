options = commandArgs(trailingOnly = TRUE)
input = options[1]
data = read.delim(input, header=T,row.names=1,sep="\t")
data = data[,(colSums(data > 0) >= 0.1)]
data = data[, grep("genus.",colnames(data))]
library(vegan)
pco1 = capscale(data ~ 0,distance = "bray")
pcs = pco1$CA$u[,1:10]
rownames(pcs) = paste0("Sample",1:nrow(pcs))
write.table(pcs, file = "10PCS.txt",sep="\t")