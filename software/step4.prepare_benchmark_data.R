options = commandArgs(trailingOnly = TRUE)
input_taxonomy = options[1]
input_annotation = paste0(options[1],".annot")
output_folder = options[2]
taxonomy_table = read.table(input_taxonomy,header=T,sep="\t")
colnames(taxonomy_table)[1] = "-"
annot_table = read.table(input_annotation,header=T,as.is = T)

taxa = c("genus.Bacteroides.id.918",
"genus.Faecalibacterium.id.2057",
"genus.Prevotella9.id.11183",
"genus.Blautia.id.1992",
"genus.Alistipes.id.968",
"genus.Subdoligranulum.id.2070",
"genus..Eubacteriumrectalegroup.id.14374",
"genus.Bifidobacterium.id.436",
"genus.Dialister.id.2183",
"genus.Ruminococcus2.id.11374")


dir.create(output_folder)

for (i in taxa){
	tax = taxonomy_table[i,]
	annot = annot_table[annot_table[,2] == i,]
	write.table(tax,file = paste0(output_folder,"/",i,".txt",sep="\t",quote=F)
	write.table(annot,file = paste0(output_folder,"/",i,".txt.annot",sep="\t",quote=F)
}
