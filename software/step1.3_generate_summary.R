options = commandArgs(trailingOnly = TRUE)
input = options[1]
output = options[2]

data = read.table(input, header=T,row.names=1,sep="\t")
summary.results = data.frame(zeros = colSums(data==0),
                             nonzeros = colSums(data>0),
                             mean.ab.reads = apply(data,2,function(x){mean(x[x>0])}),
                             mean.ab.share = apply(data,2,function(x){x = x/data[,"rootrank.Root"];mean(x[x>0])}))

write.table(summary.results,file = options[2],sep="\t")