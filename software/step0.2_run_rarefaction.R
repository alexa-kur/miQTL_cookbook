options = commandArgs(trailingOnly = TRUE)
file = options[1]
cutoff = options[2]
system(paste0("nl ",file, "|grep '>' > coord.tbl"))

coord =read.table("coord.tbl",as.is = T)
coord[1:(nrow(coord)-1),3] = coord[2:nrow(coord),1] -1 

coord[nrow(coord),3] = as.integer(sub(" .*", "", system(paste0("wc -l ", file),intern = T)))
coord[,4] = factor(sub("_.*","",coord$V2))
samples = table(coord$V4)
sample_names = names(samples)
good_samples = sample_names[samples > cutoff]
result = data.frame(V1 = as.integer(),V2 = as.character(),V3 = as.integer(),V4=as.character())

set.seed(1234)
for(i in good_samples){
  subset = coord[coord$V4 == i,]
  s = sample(1:nrow(subset),size = cutoff)
  subset2 = subset[s,]
  result = rbind(result,subset2)
}
data.new = data.frame(V1 = sub(">","",result$V2))

write.table(data.new,file = "2filter.ids",row.names = F,col.names = F,quote = F)
system("rm coord.tbl")
#system("source activate qiime1;filter_fasta.py -f LLD.fna -o filter.fa -s 2filter.ids")