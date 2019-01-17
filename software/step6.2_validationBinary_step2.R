# load libraries ----------------------------------------------------------
if(!require(foreach)){
  install.packages("foreach",repos='http://cran.us.r-project.org')
  library(foreach)
}
if(!require(getopt)){
  install.packages("getopt",repos='http://cran.us.r-project.org')
  library(getopt)
}

spec = matrix (c(
  'genotypes_folder' , 'g' , 2, 'character',
  'phenotypes_file' , 'p' , 2, 'character',
  'microbe_file' , 'm' , 2, 'character',
  'results_file' , 'r' , 2, 'character',
  'coupling_file' , 'c' , 2, 'character',
  'harmonizer' , 'h' , 2, 'character',
  'name', 'n', 2, 'character'
),byrow = T,ncol = 4)
opt = getopt(spec)

if(is.null(opt$genotypes_folder))   {opt$genotypes_folder = "./trityper_chapter4"}
if(is.null(opt$phenotypes_file) )   {opt$phenotypes_file = "./phenotypes.txt"}
if(is.null(opt$microbe_file)    )   {opt$microbe_file = "./tax_binary.txt"}
if(is.null(opt$results_file)    )   {opt$results_file = "./topBinary_Jan2019.txt"}
if(is.null(opt$coupling_file)   )   {opt$coupling_file = "./coupling_file.txt"}
if(is.null(opt$harmonizer)      )   {opt$harmonizer = "./GenotypeHarmonizer-1.4.20-SNAPSHOT"}
if(is.null(opt$name)            )   {message("BINARY VALIDATION SCRIPT MESSAGE: Please define the name of your cohort: Rscript -n COHORT");q(status = 1)}

if(!file.exists(opt$genotypes_folder)){message("genotypes folder doesn't exist!");quit(save="no")}
if(!file.exists(opt$phenotypes_file)){message("phenotypes file doesn't exist!");quit(save="no")}
if(!file.exists(opt$microbe_file)){message("microbe file doesn't exist!");quit(save="no")}
if(!file.exists(opt$results_file)){message("results file folder doesn't exist!");quit(save="no")}
if(!file.exists(opt$coupling_file)){message("coupling file folder doesn't exist!");quit(save="no")}
if(!file.exists(opt$harmonizer)){message("path to GenotypeHarmonizer folder is not correct!");quit(save="no")}


.Last = function(){
  
  if(file.exists(paste0(opt$genotypes_folder,"/SNPs.txt.backup"))){ 
    system(paste0("mv ",opt$genotypes_folder,"/SNPs.txt.backup ", opt$genotypes_folder,"/SNPs.txt"))
    }
  if(file.exists(paste0(opt$genotypes_folder,"//SNPMappings.txt.backup")))   {
    system(paste0("mv ",opt$genotypes_folder,"/SNPMappings.txt.backup ", opt$genotypes_folder,"/SNPMappings.txt"))
  }
     
}



message("all necessary files exist, trying to load SNP subset")

# Creating dosage files ---------------------------------------------------
data = read.table(opt$results_file,header=T,as.is = T,sep="\t")
data$ST.bac.name = sub(".*id","id",data$HGNCName)

full_snpnames = read.table(paste0(opt$genotypes_folder,"/SNPs.txt"),header=F,as.is = T)

if(any(duplicated(full_snpnames[,1]))){
  if(!file.exists(paste0(opt$genotypes_folder,"/SNPs.txt.backup"))) {
    system(paste0("cp ",opt$genotypes_folder,"/SNPs.txt ", opt$genotypes_folder,"/SNPs.txt.backup"))
    } else {message("backup file exists already");quit(save="no")}
  if(!file.exists(paste0(opt$genotypes_folder,"/SNPMappings.txt.backup"))) { 
    system(paste0("cp ",opt$genotypes_folder,"/SNPMappings.txt ", opt$genotypes_folder,"/SNPMappings.txt.backup"))
  } else {message("backup file exists already");quit(save="no")}
  full_snplist2update = read.table(paste0(opt$genotypes_folder,"/SNPMappings.txt"),as.is = T,sep="\t")
  
  full_snpnames[duplicated(full_snpnames[,1]),1] = paste0(full_snpnames[duplicated(full_snpnames[,1]),1],".2")
  full_snpnames[duplicated(full_snpnames[,1]),1] = paste0(full_snpnames[duplicated(full_snpnames[,1]),1],".3")
  full_snplist2update[,3] = full_snpnames[,1]
  write.table(full_snpnames, file = paste0(opt$genotypes_folder,"/SNPs.txt"),col.names=F,quote = F,row.names=F,sep="\t")
  write.table(full_snplist2update, file = paste0(opt$genotypes_folder,"/SNPMappings.txt"),col.names=F,quote = F,row.names=F,sep="\t")
}
rm(full_snpnames)
rm(full_snplist2update)

full_snplist = read.table(paste0(opt$genotypes_folder,"/SNPMappings.txt"),as.is = T,sep="\t")

full_snplist$standard_name = as.character(paste(full_snplist$V1,full_snplist$V2,sep=":"))

#full_snplist = full_snplist[!duplicated(full_snplist$standard_name),]

snp_selection = full_snplist$V3[full_snplist$standard_name %in% data$SNPName]

if (length(snp_selection) < 1){
  message("There's a problem on SNP selection stage. e-mail Alex Kurilshikov (alexa.kur@gmail.com) to solve the issue")
  quit(save = "no")
}
write.table(snp_selection,"tmp.snp.subset.txt",quote = F,row.names=F,col.names=F)

message("external call of genotypeHarmonizer")

if (file.exists(paste0(opt$harmonizer,"/GenotypeHarmonizer.jar"))){
  system(paste0("java -Xmx40G -jar ", opt$harmonizer, "/GenotypeHarmonizer.jar -i ",
                opt$genotypes_folder, 
                " -I TRITYPER -o trityper_subset -O TABLE -vf ./tmp.snp.subset.txt"))
} else {
  message ("please put GenotypeHarmonizer-1.4.20-SNAPSHOT in the current folder or define a path to its folder using -h parameter!")
  quit(save = "no")
}

message("SNP subset preparation done. Harmonizing data...")
# Loading and harmonizing all the data ------------------------------------

dosages = read.table("trityper_subset.dosages.txt",header=T,row.names=1,sep="\t",check.names = F)
genotypes = read.table("trityper_subset.genotypes.txt",header=T,row.names=1,sep = "\t",check.names =F)
coupling = read.table(opt$coupling_file,header=F,as.is = T,colClasses = "character")
phenotypes = read.table(opt$phenotypes_file,row.names=1,header=T,sep="\t")
microbe = read.table(opt$microbe_file,header=T,row.names=1,sep="\t",check.names = F) - 100

paired_samples = (coupling[,1] %in% colnames(dosages)) & (coupling[,2] %in% colnames(microbe)) & (coupling[,2] %in% rownames(phenotypes))

coupling_subset = coupling[paired_samples,]

dosages = t(dosages)[coupling_subset[,1],]
genotypes = t(genotypes)[coupling_subset[,1],]
microbe = t(microbe)[coupling_subset[,2],]
colnames(microbe) = sub(".*id","id",colnames(microbe))
phenotypes = phenotypes[coupling_subset[,2],]

rownames(full_snplist) = full_snplist$V3
#colnames(dosages) = full_snplist[colnames(dosages),]$standard_name
#colnames(genotypes) = full_snplist[colnames(genotypes),]$standard_name
index.snpnames = full_snplist[colnames(dosages),]$standard_name

message(paste0("all data layers loaded. total number of samples for analysis: ", nrow(dosages)))
# running_model and report results ----------------------------------------

foreach(i = 1:nrow(data),.combine = rbind) %do% {

  if ((data$ST.bac.name[i] %in% colnames(microbe)) & (data$SNPName[i] %in% index.snpnames) ) {
    bac = microbe[,data$ST.bac.name[i]]
    dos = dosages[,index.snpnames == data$SNPName[i],drop = F]
    gen = genotypes[,index.snpnames == data$SNPName[i],drop = F]
   result = foreach(j = 1:ncol(dos),.combine = rbind) %do%{
      glm1 = glm(bac ~ dos[,j] + .,data = phenotypes,family = "binomial")
      summary1 = summary(glm1)
      data.frame(bac = data$ST.bac.name[i],
                 snp = data$SNPName[i],
                 ref.allele = gen[which.min(dos[,j]),j],
                 alt.allele = gen[which.max(dos[,j]),j],
                 beta = summary1$coef[2,1],
                 sd = summary1$coef[2,2],
                 z = summary1$coef[2,3],
                 p = summary1$coef[2,4],
                 nsamples  = summary1$df.null + 1
      )
    }
    
  } else {
    result = data.frame(bac = data$ST.bac.name[i], snp = data$SNPName[i], ref.allele = NA, alt.allele = NA, beta = NA,sd = NA, z= NA, p = NA, nsamples = NA)
  }
  result
}-> report

write.table(report, file = paste0(opt$name,"_binaryResults_Jan2019_updated.txt"),sep="\t",quote = F,row.names=F)
system(paste0("gzip ", opt$name,"_binaryResults_Jan2019_updated.txt"))

#system(paste0("mv ",opt$genotypes_folder,"/SNPs.txt.backup ", opt$genotypes_folder,"/SNPs.txt"))
#system(paste0("mv ",opt$genotypes_folder,"/SNPMappings.txt.backup ", opt$genotypes_folder,"/SNPMappings.txt"))


system("rm trityper_subset.dosages.txt trityper_subset.genotypes.txt tmp.snp.subset.txt")
message(paste0("Finished succesfully. Please send the file ",opt$name,"_binaryResults_updated.txt.gz to Alex Kurilshikov: alexa.kur@gmail.com"))
quit(save="no",runLast = TRUE)





