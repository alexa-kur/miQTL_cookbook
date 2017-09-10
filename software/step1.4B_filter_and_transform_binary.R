options = commandArgs(trailingOnly = TRUE)
input = options[1]
output = options[2]
cutoff = as.numeric(options[3])
data = read.table(input, header=T,row.names=1,sep="\t")
data = data[,(colSums(data > 0) >= cutoff*nrow(data))]
data = data[, grep("NOTAX",colnames(data),invert = T)]
data = data[,grep("domain",colnames(data),invert = T)]
data = data[,grep("rootrank.Root",colnames(data),invert = T)]
data[data == 0] = NA
data = log(data)
write.table(data,file = options[2],sep="\t",quote = F)
