library(foreach)

#usage: Rscript model_prototype.R microbiome.txt genotypes.txt phenotypes.txt
options = commandArgs(trailingOnly = TRUE)
microbiome = read.table(options[1],header=T,row.names=1,sep="\t")
genotypes = t(read.table(options[2],header=T,row.names=1,sep="\t"))
phenotypes = read.table(options[3],header=T,row.names=1,sep="\t")

common_names = intersect(rownames(microbiome),
                         intersect(rownames(genotypes),
                                   rownames(phenotypes)))

message(paste0("Number of samples with all data types available: ",length(common_names)))

microbiome = microbiome[common_names,]
genotypes = genotypes[common_names,]
phenotypes = phenotypes[common_names,]

foreach (i = 1:ncol(microbiome),.combine = rbind)%do%{
  print(i)
  tax_vector = microbiome[,i]
  tax_vector_bin = as.integer(!is.na(microbiome[,i]))
  binary_zero = glm(tax_vector_bin ~ . ,data = phenotypes,family = "binomial")
  quant_zero = lm(tax_vector ~ .,data = phenotypes)
  
  foreach(j = 1:ncol(genotypes),.combine = rbind)%do%{
    if( j %% 100 == 0) print(j) 
    snp = genotypes[,j]
    binary = glm(tax_vector_bin ~ . + snp,data = phenotypes,family = "binomial")
    quant = lm(tax_vector ~ . + snp,data = phenotypes)
    anova.bin = anova(binary_zero,binary,test = "Chisq")
    anova.quant = anova(quant_zero,quant)
    data.frame(trait = colnames(microbiome)[i],
               snp = colnames(genotypes)[j],
               coef.bin = binary$coef[length(binary$coef)],
               p.bin = anova.bin[2,5],
               coef.quant = quant$coef[length(quant$coef)],
               p.quant=anova.quant[2,6])
  }
}



foreach (i = 1:ncol(microbiome),.combine = rbind)%do%{
  print(i)
  tax_vector = microbiome[,i]
  tax_vector_bin = as.integer(!is.na(microbiome[,i]))
  binary_zero = glm(tax_vector_bin ~ . ,data = phenotypes,family = "binomial")
  quant_zero = lm(tax_vector ~ .,data = phenotypes)
  
  
  result = apply(genotypes[,1:1000],2,function(snp){
    binary = glm(tax_vector_bin ~ . + snp,data = phenotypes,family = "binomial")
    quant = lm(tax_vector ~ . + snp,data = phenotypes)
    anova.bin = anova(binary_zero,binary,test = "Chisq")
    anova.quant = anova(quant_zero,quant)
    return(c(
               coef.bin = binary$coef[length(binary$coef)],
               p.bin = anova.bin[2,5],
               coef.quant = quant$coef[length(quant$coef)],
               p.quant=anova.quant[2,6]))
  })
  t(result)
  #result = as.data.frame(t(result))
  #data.frame(bac = colnames(microbiome)[i],snp = colnames(genotypes)[1:10000],result)
}-> result



