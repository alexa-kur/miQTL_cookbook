options = commandArgs(trailingOnly = TRUE)
taxonomy_file = options[1]
output_file = options[2]

tax = read.table(taxonomy_file,sep="\t",header=T,as.is = T)

annot = data.frame(Platform="RDP",HT12v3.ArrayAddress=tax[,1],Gene=tax[,1],Chr=4,ChrStart=1000,ChrEnd=1000)
colnames(annot)[2]="HT12v3-ArrayAddress"

write.table(annot,file = output_file,sep="\t",quote = F,row.names=F)