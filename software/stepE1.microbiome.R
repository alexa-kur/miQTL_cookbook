#Multiple microbiome stats using the microbiome package (https://microbiome.github.io/microbiome/)

#Load libraries
library("phyloseq"); packageVersion("phyloseq")
require(microbiome)

if(!require(BiocInstaller)){
  install.packages("foreach",repos='http://cran.us.r-project.org')
  library(BiocInstaller)
}

if(!require(microbiome)){
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("microbiome")
  library(microbiome)
}
if(!require(phyloseq)){
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("phyloseq")
  library(phyloseq)
}

if(!require(getopt)){
  install.packages("getopt",repos='http://cran.us.r-project.org')
  library(getopt)
}
spec = matrix (c(
  'input' , 'i' , 2, 'character',
  'output','o',2,'character'
),byrow = T,ncol = 4)
opt = getopt(spec)

if(is.null(opt$input))   {opt$input = "./taxonomy_table.txt"}
if(is.null(opt$input))   {opt$output = "./microbiomeScript_summary"}

if(!file.exists(opt$input)) {message("taxonomy_table not provided or not exist!");q(status = 1)}
if(dir.exists(opt$output)) {message("output folder exists! please delete before running the script");q(status = 1)}
if(!dir.exists(opt$output)) {dir.create(opt$output)}

#Read taxonomical table (only genera needed)
table <- read.table(opt$input, header=TRUE,row.names=1,sep="\t")
table = table[,grep("^genus[.]",colnames(table))]
table = table[sample(1:nrow(table)),]
rownames(table) = paste0("Sample",1:nrow(table))
table = as.data.frame(t(table))
#Convert to pseq and ensure taxa are rows
genustable <- otu_table(table, taxa_are_rows = TRUE)

genustable[1:3,1:5]
#OTU Table:          [3 taxa and 5 samples]
#taxa are rows
#Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
#G_Abiotrophia             0      0       0       0       0
#G_Acetanaerobacterium     0      0       0       0       2
#G_Acetatifactor           0      0       0       0       0

#Compute various diversities (simpson, shannon,...)
gdiversities <- diversities(genustable, index = "all")
write.csv(gdiversities, file = paste0(opt$output,"/gdiversities.csv"))

#Estimate richness at various detection thresholds
grichness <- richness(genustable)
write.csv(grichness, file = paste0(opt$output,"/grichness.csv"))

#Obtain dominance indices for most abundant taxa in each sample
gdominance <- dominance(genustable, index='all')
write.csv(gdominance, file = paste0(opt$output,"/gdominance.csv"))

#Quantify rare and low abundant taxa indices
grarity <- rarity(genustable, index = "all")
write.csv(grarity, file = paste0(opt$output,"/grarity.csv"))

#Compute coverage (number of taxa covering )
gcoverage <- coverage(genustable, threshold = 0.9)
write.csv(gcoverage, file = paste0(opt$output,"/gcoverage.csv"))

#Get core abundances
gcoreabundance <- core_abundance(genustable, detection = .1/100, prevalence = .95)
write.csv(gcoreabundance, file = paste0(opt$output,"/gcoreabundance.csv"))

#Inequality
ggini <- inequality(genustable)
write.csv(ggini, file = paste0(opt$output,"/ggini.csv"))

#ComputeEvenness
geven <- evenness(genustable, "all")
write.csv(geven, file = paste0(opt$output,"/geven.csv"))

