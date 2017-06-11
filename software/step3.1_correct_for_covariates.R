options = commandArgs(trailingOnly = TRUE)
taxonomy_file = options[1]
phenotype_file = options[2]
output_file = options[3]
coupling_file = options[4]

tax = read.table(taxonomy_file,header=T,row.names=1,sep="\t")
phenos = read.table(phenotype_file,header = T,row.names=1,sep="\t")

if (exists("coupling_file")){
  coupling = read.table(coupling_file,as.is = T)
  has_both = (coupling[,1] %in% rownames(phenos)) & (coupling[,2] %in% rownames(tax))
  coupling= coupling[has_both,]
  tax = tax[coupling[,2],]
  phenos = phenos[coupling[,1],]
  rownames(tax) = rownames(phenos)
} else {
  tax = tax[intersect(rownames(tax),rownames(phenos)),]
  phenos = phenos[rownames(tax),]
}

corrected_data = apply(tax,2,function(x){
  x.subset = x[!is.na(x)]
  phenos.subset = phenos[!is.na(x),]
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

write.table(corrected_data2, file = output_file,sep="\t",row.names = F,quote = F)
write.table(annot,file = paste0(output_file,".annot"),sep="\t",row.names=F,quote = F)