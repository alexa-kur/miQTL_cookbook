options = commandArgs(trailingOnly = TRUE)
taxonomy_file = options[1]
phenotype_file = options[2]
coupling_file = options[3]

if(!require(vegan)){
    install.packages("vegan")
    library(vegan)
}

tax = read.table(taxonomy_file,header=T,row.names=1,sep="\t")
phenos = read.table(phenotype_file,header = T,row.names=1,sep="\t")
phenos = phenos[!apply(phenos,1,function(x){any(is.na(x))}),,drop=F]


if (exists("coupling_file")&file.exists(coupling_file)){
  coupling = read.table(coupling_file,colClasses = "character")
  has_both = (coupling[,1] %in% rownames(phenos)) & (coupling[,2] %in% rownames(tax))
  coupling= coupling[has_both,]
  tax = tax[coupling[,2],]
  phenos = phenos[coupling[,1],,drop = FALSE]
  rownames(tax) = rownames(phenos)
  tax = tax[rownames(tax)[grep("[.][0-9]+$",rownames(tax),invert=T)],]
  phenos = phenos[rownames(tax),,drop = FALSE]
} else {
  tax = tax[intersect(rownames(tax),rownames(phenos)),]
  phenos = phenos[rownames(tax),,drop = FALSE]
}

tax = tax[,grep("^genus[.]",colnames(tax))]


shannon = diversity(tax,index = "shannon")
simpson = diversity(tax,index = "simpson")
invsimpson = diversity(tax,index = "invsimpson")
all_div = data.frame(row.names = rownames(tax),shannon = shannon, simpson = simpson,invsimpson = invsimpson)


corrected_data = apply(all_div,2,function(x){
  x.subset = x[!is.na(x)]
  phenos.subset = phenos[!is.na(x),,drop = FALSE]
  phenos.subset.matrix = data.matrix(phenos.subset)
  if(ncol(phenos.subset)==ncol(phenos.subset.matrix)){
    phenos.subset = phenos.subset[,apply(phenos.subset.matrix,2,sd) !=0,drop = FALSE]
  }
  
  x.resid = resid(lm(x.subset ~ .,data = phenos.subset))
  x[!is.na(x)] = x.resid+100
  x[is.na(x)] = 0
  return(x)
})

corrected_data = as.data.frame(t(corrected_data))
corrected_data2 = data.frame(rownames(corrected_data),corrected_data)
colnames(corrected_data2)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data2),
                   Gene = rownames(corrected_data2),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"


write.table(corrected_data2, file = "alpha_div.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "alpha_div.txt.annot",sep="\t",row.names=F,quote = F)

