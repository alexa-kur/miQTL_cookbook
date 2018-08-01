options = commandArgs(trailingOnly = TRUE)
input = options[1]
phenos_input = options[2]
coupling_file = options[3]

data = read.delim(input, header=T,row.names=1,sep="\t")
data = data[,(colSums(data > 0) >= 0.1)]
data = data[, grep("genus.",colnames(data))]

if(!require(vegan)){
    install.packages("vegan")
    library(vegan)
}

pco1 = capscale(data ~ 0,distance = "bray")
pcs = pco1$CA$u[,1:10]
rownames(pcs) = paste0("Sample",1:nrow(pcs))
write.table(pcs, file = "10PCS.txt",sep="\t")
phenos = read.delim(phenos_input, header=T,row.names=1,sep="\t")

if (exists("coupling_file")&file.exists(coupling_file)){
  coupling = read.table(coupling_file,as.is = T)
  has_both = (coupling[,1] %in% rownames(phenos)) & (coupling[,2] %in% rownames(data))
  coupling= coupling[has_both,]
  tax = data[coupling[,2],]
  phenos = phenos[coupling[,1],,drop = FALSE]
  rownames(tax) = rownames(phenos)
  tax = tax[rownames(tax)[grep("[.][0-9]+$",rownames(tax),invert=T)],]
  phenos = phenos[rownames(tax),,drop = FALSE]
} else {
  tax = data[intersect(rownames(data),rownames(phenos)),]
  phenos = phenos[rownames(tax),,drop = FALSE]
}

dist.mat = vegdist(tax,method = "bray")
results = list()
for (i in 1:ncol(phenos)){
  ad1 = adonis(dist.mat ~ phenos[,i],permutations = 1000)
  results[[i]] = ad1$aov.tab[1,]
}

ad_results = do.call(rbind,results)
rownames(ad_results) = colnames(phenos)
write.table(ad_results,file = "adonis_covariates.txt",sep="\t")
