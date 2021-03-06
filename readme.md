---
title: "MiBioGen miQTL pipeline"
author: "Alexander Kurilshikov"
date: "Jan 31, 2018"
output: html_document
---

## Contacts

If you have any questions and suggestions, don't hesitate to write us: Alexander Kurilshikov (alexa.kur@gmail.com), Sasha Zhernakova (sasha.zhernakova@gmail.com)

# Overall description

This is the cookbook for performing the GWAS analysis of microbial abundance based on analysis of 16S rRNA sequencing dataset. It includes 4 major steps:

1. [Processing of 16S data](#chapter-1-16s-data-processing)
2. [Processing of SNP microarray data](#chapter-2-genotype-imputation)
3. [Performing the association study benchmark](#chapter-3-perform-genome-wide-association-study-benchmark)
4. [Performing the association study for all taxa](#chapter-4-perform-genome-wide-association-study-for-all-taxa)
5. [Run GWAS on alpha diversity](#chapter-5-run-gwas-on-alpha-diversity)
6. [Run binary trait validation](#chapter-6-run-binary-trait-validation)
7. [Acquire account and upload data](#chapter-7-acquire-guest-account-and-upload-data) 
8. Making the meta-analysis

Steps from 1 to 3 would be performed in-house by every participating group. Step 4 will be performed in UMCG (Groningen). Mostly, you can just copy the code strings from this cookbook and run them, but sometimes it's not the case (yet, but it will be after beta-testing).

#### Genome-Wide Association Study itself will be performed according to this design:

1. Cutoffs and transformations
    + Taxonomies:
        + Abundance cutoff: presence in 10% of the samples
        + Log (base **e**) trasformation on the counts
    + SNPs
        + MAF > 5%
        + Imputation quality > 0.4
        + SNP call rate > 0.95
        + Genotypes represented in dosages 
2. Models used
    + Taxonomy absence/presence as binary trait: logistic regression with Chisquared-based p-value estimation
    + For non-zero samples: linear regression model on log-transformed counts with Fisher test-based p-value estimation
3. Meta-analysis
    + will be performed separately, for binary and quantitative models

# Chapter 1. 16S data processing

For 16S analysis, RDP Classifier will be used instead of OTU picking, since it has shown more consistent results between different domains. According to the pipeline, Genome-Wide Association Study will be performed for the following taxonomic levels:

1. Taxonomic levels:
    + Classes
    + Orders
    + Phyla
    + Families
    + Genera
2. Cutoffs:
    + Taxonomy should be presented in more than 10% of the samples in cohort
3. Transformations:
    + Taxonomies counts should be log-transformed (on the base ***e***)

## 1.1. 16S data requirements

1. Sequence quality-based read filtering should be performed before rarefaction.
2. Reads should be randomly rarefied to 10,000 reads per sample before OTU picking (See [Appendix2](#Appendix2) )
3. Sequences from all samples should be merged into one file.
4. For every sequence, FASTA header should follow this format:

```
>[SampleID]_[SequenceID] [OptionalInfo]
```

Example of valid fasta record:

```
>G36899_G36899.A2137
TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGACTGG
CAAGTCTGATGTGAAAGGCGGGGGCTCAACCCCTGGACTGCATTGGAAACTGTTAGTCTT
GAGTGCCGGAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAG
GAGCACCAGTGGCGAAGGCGGCTTACTGGACGGTAACTGACGTTGAGGCTCGAAAGCGTG
GGGAGCAAACAGG
```

For this record, SampleID is **G36899**, and SequenceID is **G36899.A2137.451**.

## 1.2. Software and database installation

These software and databases are necessary:

1. Java
2. R
3. RDP Classifier 2.12. Can be downloaded [here](https://sourceforge.net/projects/rdp-classifier/)

You also need reference database for RDP Classifier and some additional scripts. They are included in this GitHub page, and located in the folders [database](https://github.com/alexa-kur/miQTL_cookbook/tree/master/database) and [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software). The list of things you need to get here includes:

4. SILVA 128 release database prepared for RDP Classifier. available on this Github in the folder named [database](https://github.com/alexa-kur/miQTL_cookbook/tree/master/database). It's a spanned zip archive, so you first need to unzip it. It's also available on the consortium [Dropbox](https://www.dropbox.com/home/Microbiome-QTL_Charge): located in the folder [CookBook_CHARGE_miQTL/database](https://www.dropbox.com/home/Microbiome-QTL_Charge/CookBook_CHARGE_miQTL/database) and named *silva128_rdpDB.zip*. 
5. RDP output parsing script. Located in the folder [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software) and named **step1.2_rdp_process.sh**.
6. Script for generating summary statistics **step1.3_generate_summary.R**. Located in the folder [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software)
7. Filtering and transformation script **step1.4_filter_and_transform.R**. Located in the folder [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software)

Before usage, you can put the scripts and 16S database folter to your project folder. 

## 1.3. Running data processing. 

### Step 1.1. Taxonomy binning 

Go to your project folder. Replace SEQUENCES.FASTA to the filename of your 16S sequences and run:

```
java -Xmx10G -jar ./rdp_classifier_2.12/dist/classifier.jar -t ./SILVA128_rdpDB/rRNAClassifier.properties -o results.out SEQUENCES.FASTA
````

File **results.out** will be generated. Please note that this step is time consuming and can take up to several days. 


### Step 1.2. Process mapping results

```
bash step1.2_rdp_process.sh results.out 0.8
````

It will generate file **taxonomy_table.txt**, tab-separated table which contains per-sample counts of taxa in your dataset. 


### Step 1.3. Generate summary statistics

Please run this code and send us the results file **COHORT_NAME_summary_16s.txt** (replace COHORT_NAME to the real name of your cohort)

```
Rscript step1.3_generate_summary.R taxonomy_table.txt COHORT_NAME_summary_16s.txt
```

### Step 1.4. Abundance filtering and transformation

The script *step1.4_filter_and_transform.R* removes taxa presented in less than 10% of samples and performes log transformation 
```
Rscript step1.4_filter_and_transform.R taxonomy_table.txt tax_filtered_logTrans.txt 0.1
```

**tax_filtered_logTrans.txt** fill be generated. This file will be used in [miQTL mapping](#chapter3-peform-genome-wide-association-study) .  

# Chapter 2. Genotype imputation

To simplify the analysis, we should have genotypes harmonized between cohorts. We propose to use one imputation server: [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html). 
Carolina Medina-Gomez was very pleased to prepase very nice instruction for imputation, please find here, in the file called [Imputation_protocol_description.docx](Imputation_protocol_description.docx). In addition to that, you can also check the documentation on HRC [server help page](https://imputationserver.sph.umich.edu/index.html#!pages/help).

# Chapter 3. Perform Genome-Wide Association Study Benchmark

Compare to the whole microbiome analysis prepared on the next step, here we'll only focus on 10 most abundant bacterial genera. This will allow us to explore cohort-specific stange

In this section, we will perform the following steps:

1. [Correct microbiome data for covariates](#step-31-correct-microbiome-for-covariates-and-create-annotation-file)
2. [Process HRC VCF files](#step-32-process-hrc-vcf-files)
3. [Prepare coupling file](#step-33-prepare-coupling-file)
4. [Rename cohort in template](#step-34-rename-dataset-in-template)
5. [Run miQTL mapping](#step-34-run-miqtl-mapping)

All necessary scripts are provided in [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software) folder. Thus, it includes:

1. Script **step3.1_correct_for_covariates.R**
2. **GenotypeHarmonizer**
3. **eQTL-pipeline**

## Step 3.1. Correct microbiome for covariates and create annotation file

**Important note 1: Please, re-do Step 1.4 with the most updated script!**

**Important note 2: Please, for initial run perform the correction only for age and gender!**

First, you need to generate filtered and log-transformed taxonomy table *tax_filtered_logTrans.txt*. Please check [Step 1.4](#step-14-abundance-filtering-and-transformation) to explore how to do it. Now you need to create the flat text table with phenotypes you want to correct for. It should be tab-delimited flat text file which suits for loading in R. Phenotypes can be integer or numeric variables or factors. The first column should contain Sample ID's. The typical example is provided below: 

```
SampleID  Sex   Age
Sample1   0       23
Sample2   1       74
Sample3   0       58
Sample4   1       76
```

This also works (column name for SampleID's can be skipped)

```
Sex   Age
Sample1   0       23
Sample2   1       74
Sample3   0       58
Sample4   1       76
```

If your sample names are **different** between taxonomy table and phenotype file, please add **linkage file**. This is an example of linkage file:

```
PhenotypeSample1  StoolSample1
PhenotypeSample2  StoolSample2
```

Phenotype ID comes first, followed by stool sample ID (again, it's tab-delimited file). 

Now run (assuming your phenotype file is called **phenotypes.txt**)

```
#If your sample ID's are matched
Rscript step3.1_correct_for_covariates.R tax_filtered_logTrans.txt phenotypes.txt tax_corrected.txt

#If your sample ID's are different between phenotypes and 16s libraries
Rscript step3.1_correct_for_covariates.R tax_filtered_logTrans.txt phenotypes.txt tax_corrected.txt linkage.txt
```

The files **tax_corrected.txt** and **tax_corrected.txt.annot** will be generated. 

## Step 3.1B. Prepare benchmarking data
Next, we need to generate individual files for each of 10 bacterial genera for benchmarking purpose.

```
Rscript step3.1B_prepare_benchmark_data.R tax_corrected.txt 
```

it will create a folder **taxa_benchmark_selection** which contain all necessary files for running GWAS on 10 bacteria

## Step 3.2. Process HRC VCF files

First, you need to install **GenotypeHarmonizer** software. simply untar the archive [GenotypeHarmonizer-1.4.20-dist.tar.gz](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/GenotypeHarmonizer-1.4.20-dist.tar.gz) from [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software) folder:

```
tar -zxf GenotypeHarmonizer-1.4.20-dist.tar.gz
``` 

Please check the structure of the genotype folder you downloaded from HRC. It should be per-chromosome gzipped VCF files:

```
chr10.dose.vcf.gz
chr10.dose.vcf.gz.md5
chr10.dose.vcf.gz.tbi
chr10.dose.vcf.gz.tbi.md5
chr10.info.gz
chr10.info.gz.md5
chr11.dose.vcf.gz
chr11.dose.vcf.gz.md5
chr11.dose.vcf.gz.tbi
chr11.dose.vcf.gz.tbi.md5
...
```

Keep only autosomal data in the folder (assuming your HRC folder is named **HRC_VCF**):

```
mkdir autosomal_genotypes
cp HRC_VCF/chr[0-9]* autosomal_genotypes
```

Next command will perform the following steps:
1. MAF filtering (*-mf* flag, we decided to use 0.05 threshold)
2. Imputation QC filtering (*-ip* flag, 0.4)
3. Remove ambigous SNPs (*-asf* flag)
4. Translate genotypes from VCF to TRITYPER format

To apply these procedures, run:
```
java -Xmx40G -jar ./GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i ./autosomal_genotypes -I VCFFOLDER -o genotypes_trityper -O TRITYPER -mf 0.05 -ip 0.4 -asf  
```

The folder **genotypes_trityper** will be generated.

For further information on *GenotypeHarmonizer*, please check its [GitHub](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) or just contact me (alexa.kur@gmail.com).


## Step 3.3. Prepare coupling file

This file is necessary to run the pipeline. It should contain linkages between genotype and microbiome sample names. It's necessary even if you have the same names (in this case, just make two exact columns). It should be tab-delimited flat text file with two columns:

```
GenotypeSample1  StoolSample1
GenotypeSample2  StoolSample2
```


## Step 3.4. Rename dataset in templates

Our eQTL pipeline requires all the configuraiton written in XML format. The [template benchmark folder contains several XMLs](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/benchmark_templates) can be downloaded from [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software), for one genus each.

If you used standard file and folder names mentioned in this cookbook, the only change you need to do in each template is the dataset name. that's easy one-step procedure. For example, assuming your cohort is named *MBRUN*, run

```
for i in `ls miQTL_cookbook/software/benchmark_templates` ; do j=${i/template_/};cat miQTL_cookbook/software/benchmark_templates/${i} |perl -pe 's/COHORTNAME/MBRUN/' > ${j} ; done
```

it will generate all 10 XML benchmarks for further analysis. 

## Step 3.5. Run miQTL mapping per each of 10 genera

First you need to unzip eQTL pipeline binaries: [eqtl-mapping-pipeline-1.4nZ-dist.zip](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/eqtl-mapping-pipeline-1.4nZ-dist.zip) from [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software). Simply run:

```
unzip eqtl-mapping-pipeline-1.4nZ-dist.zip
```
These are the files that should be located in project folder:

1. **tax_corrected.txt**  // The microbiome file which contains filtered, transformed and covariate-corrected microbial taxa. See [Chapter 1](#chapter-1-16s-data-processing) and [Step 3.1.](#step-31-correct-microbiome-for-covariates-and-create-annotation-file)
2. **tax_corrected.txt.annot** //Annotation file for taxonomies. See [Step 3.1.](#step-31-correct-microbiome-for-covariates-and-create-annotation-file)
3. **genotypes_trityper** // The folder which contains MAF and QC filtered autosomal genotypes in TriTyper format. See [Step 3.2](#step-32-process-hrc-vcf-files)
4. **coupling_file.txt** // Coupling file. See [Step 3.3.](#step-33-prepare-coupling-file)
5. **eqtl-mapping-pipeline-1.4nZ** // The folder with eQTL pipeline binaries
6. 10 files called **benchmarkN.xml** // XML settings files. See [Step 3.4.](#step-34-rename-dataset-in-templates)

If you have all of these ready, for the first behcnmark (**benchmark0.xml**) template you can just run this from a *screen* session:

```
java -XX:ParallelGCThreads=5 -Xmx30G -jar eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings benchmark0.xml 
```

Repeat this procedure for all 10 files. You will see your output inf the folder named **MBRUN_output**. In each subfolder you will find all QTL effects in files **eQTLs.txt.gz**. 

In addition to playing with these results we want you to upload it to UMCG servers. Please go to [Chapter 7](#chapter-7-acquire-guest-account-and-upload-data) to find how to get access and upload the data. Please note that we need all the files within output folder.

## Step 3.6. PCoA plot on genera level

We need to explore what kind of technical covarites significantly affect your microbiome composition. To do so, we want to perform PCoA analysis and PERMANOVA on technical covariates as traits. If your cohort was standardized and assumed no to contain batch effects, please use phenotype file that contains covariates you are going to correct for (by default, it's age and gender)

Please first check that you have *vegan* package is installed in your R. To do so, open your R interactive shell and run:

```
bash> R
R> install.packages("vegan")
R> quit
```

Next, run simple script. The usage is quite similar to step **Step 3.1**. I may or may not provide additional linkage file. The format of phenotype and linkage file should be similar to **Step 3.1**

```
#If your sample ID's are different between microbial and covariate tables
Rscript step3.5_pcoa.R taxonomy_table.txt technical_covariates.txt coupling_file.txt
#If your sample ID's are the same
Rscript step3.5_pcoa.R taxonomy_table.txt technical_covariates.txt
```

This script will generate 2 files: **10PCS.txt** and **adonis_covariates.txt**. Please send it Alex Kur *alexa.kur@gmail.com*

# Chapter 4. Perform Genome-Wide Association Study for all taxa

This Chapter describes the procedure of running GWAS on the whole list of microbiome traits. In general, the procedure is very similar to benchmark run described in [Chapter 3](#chapter-3-perform-genome-wide-association-study-benchmark). However, it contains many changes compared to Benchmark Run, so we recommend to start from the output of [Chapter 1](#chapter-1-16s-data-processing) and [Chapter 2](#chapter-2-genotype-imputation), and **avoid** using any data generated in [Benchmark run, Chapter 3](#chapter-3-perform-genome-wide-association-study-benchmark). 

To go through Chapter 4, please check that you have **all necessary files and folders**. 

1. **tax_filtered_logTrans.txt**. This file contains transformed, abundance-filtered bacterial taxa
2. **phenotypes.txt** file, the flat text tab-delimited file that contains information on covariates (see [description below](#step-41-correct-microbiome-for-covariates-create-annotation-files)). The list of necessary covariates includes:
   * **Sex**  
   * **Age**
   * **3 first genetic PCs** for monoethnic cohort OR **10 genetic PCs** for multiethnic or admixed cohort
   * **Technical microbiome covariates**, if any of these available, such as stool storage/delivery time, sequencing batches, etc. 

3. **HRC imputed genotypes** in default HRC output format, that comprises per-chromosome VCF and INFO files. See more detailed description below. 

It's also important to remove some samples. Please [Step 4.3](#step-43-prepare-coupling-file) to get the instructions: 
* **Remove ethnic outliers**, if you have any in your monoethnic cohort
* **Keep only one individual from each twin pair**, if you have cohort of twins


This chapter comprises the following steps
1. [Step 4.1. Correct microbiome for covariates, create annotation files](#step-41-correct-microbiome-for-covariate-create-annotation-files)
2. [Step 4.2. Filter genotypes](#step-42-filter-genotypes)
3. [Step 4.3. Prepare coupling file](#step-43-prepare-coupling-file)
4. [Step 4.4. Generate XML setting files for analyses](#step-44-generate-xml-setting-files-for-quantitative-and-binary-analyses)
5. [Step 4.5. Run full miQTL mapping](#step-45-run-full-miqtl-mapping)

## Step 4.1. Correct microbiome for covariates, create annotation files

The initial files for this step are  **tax_filtered_logTrans.txt** and **phenotypes.txt** . If your phenotype samples IDs and microbiome sample IDs are different you can simply add **linkage file**. 

This is an example of correct **phenotypes.txt** file. It should be *flat text*, *tab delimited*, with *header* and *sample names in rows*. If you use categorical variables (such as batches), please use *character strings* to encode them instead of *integers*. For example, if you will encode Batch1, Batch2 and Batch3 as 1, 2 and 3, the script will recognize it as numeric variable and will do the correction wrong.   

```
SampleID  Sex   Age   PC1   PC2   PC3   Batch
Sample1   0       23   1.2  9.2   0.1   b1
Sample2   1       74   3.1  1.1   3.4   b2
Sample3   0       58   2.5  -0.3  0.1   b1
Sample4   1       76   0.1  5.3   4.3   b3
```

This also works (column name for SampleID's can be skipped)

```
Sex   Age   PC1   PC2   PC3   Batch
Sample1   0       23   1.2  9.2   0.1   b1   
Sample2   1       74   3.1  1.1   3.4   b2
Sample3   0       58   2.5  -0.3  0.1   b1
Sample4   1       76   0.1  5.3   4.3   b3
```

If your sample names are **different** between taxonomy table and phenotype file, please add **linkage file**. This is an example of linkage file:

```
PhenotypeSample1  StoolSample1
PhenotypeSample2  StoolSample2
```

Phenotype ID comes first, followed by stool sample ID (again, it's tab-delimited file). 

Now run (assuming your phenotype file is called **phenotypes.txt**)

```
#If your sample ID's are matched
Rscript step4.1_covariate_correction.R tax_filtered_logTrans.txt phenotypes.txt

#If your sample ID's are different between phenotypes and 16s libraries
Rscript step4.1_covariate_correction.R tax_filtered_logTrans.txt phenotypes.txt linkage.txt
```

There will be 4 files generated:

1. **tax_numeric.txt**
2. **tax_numeric.txt.annot**
3. **tax_binary.txt**
4. **tax_binary.txt.annot**

**Important note: the samples with NA values in one or more phenotypes will be excluded from the output. If you still want to include them, please impute missing values by median or mean**

## Step 4.2. Filter genotypes

**Important note: this procedure seems very similar to the one in Benchmark run. However, for this run we will use different cutoffs, so please don't skip this step!**

First, you need to install **GenotypeHarmonizer** software. simply untar the archive [GenotypeHarmonizer-1.4.20-dist.tar.gz](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/GenotypeHarmonizer-1.4.20-dist.tar.gz) from [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software) folder:


```
tar -zxf GenotypeHarmonizer-1.4.20-dist.tar.gz
``` 

Please check the structure of the genotype folder you downloaded from HRC. It should be per-chromosome gzipped VCF files:

```
chr10.dose.vcf.gz
chr10.dose.vcf.gz.md5
chr10.dose.vcf.gz.tbi
chr10.dose.vcf.gz.tbi.md5
chr10.info.gz
chr10.info.gz.md5
chr11.dose.vcf.gz
chr11.dose.vcf.gz.md5
chr11.dose.vcf.gz.tbi
chr11.dose.vcf.gz.tbi.md5
...
```

Keep only autosomal data in the folder (assuming your HRC folder is named **HRC_VCF**):

```
mkdir autosomal_genotypes
cp HRC_VCF/chr[0-9]* autosomal_genotypes
```

Next command will perform the following steps:
1. MAF filtering (*-mf* flag, 0.05 threshold)
2. Imputation QC filtering (*-ip* flag, 0.4)
3. Call rate filtering (*-cf* flag, 0.95)
4. Translate genotypes from VCF to TRITYPER format

To apply these procedures, run:
```
java -Xmx40G -jar ./GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i ./autosomal_genotypes -I VCFFOLDER -o trityper_chapter4 -O TRITYPER -mf 0.05 -ip 0.4 -cf 0.95 
```

The folder **trityper_chapter4** will be generated.

For further information on *GenotypeHarmonizer*, please check its [GitHub](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) or just contact me (alexa.kur@gmail.com).

## Step 4.3. Prepare coupling file

**Note: if you have participated in benchmark run, you probably have this file already.**

This file is necessary to run the pipeline and should be named **coupling_file.txt**. It should contain linkages between genotype and microbiome sample names. It's necessary even if you have the same names (in this case, just make two exact columns). It should be the tab-delimited flat text file with two columns:

```
GenotypeSample1  StoolSample1
GenotypeSample2  StoolSample2
```

**If you need to exclude some samples from analysis**, the easiest way to do so is just to remove corresponding lines from this **coupling file**.

## Step 4.4. Generate XML setting files for quantitative and binary analyses

To run whole-genome miQTL mapping for all bacterial traits, we need to generate two XML setting files, one for quantitative analysis and one for bacterial presence/absence analysis. The [template XML for quantitative analysis](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/template.xml) and [template XML for binary analysis](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/template_binary.xml) can be downloaded from [software folder of this Cookbook](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software).

If you used standard file and folder names mentioned in this cookbook, the only change you need to do in templates is the dataset name. That's easy one-step procedure for each template. For example, assuming your cohort is named *MBRUN*, run

```
cat miQTL_cookbook/software/template.xml|perl -pe 's/COHORTNAME/MBRUN/' > ./settings.xml
cat miQTL_cookbook/software/template_binary.xml|perl -pe 's/COHORTNAME/MBRUN/' > ./settings_binary.xml
```

The files **settings.xml** and **settings_binary.xml** will be used on the next stage.

## Step 4.5. Run full miQTL mapping

First, you need to unzip eQTL pipeline binaries: [eqtl-mapping-pipeline-1.4nZ-dist.zip](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software/eqtl-mapping-pipeline-1.4nZ-dist.zip) from [software](https://github.com/alexa-kur/miQTL_cookbook/tree/master/software). Simply run:

```
unzip eqtl-mapping-pipeline-1.4nZ-dist.zip
```
These files that should be located in your project folder:

1. **tax_numeric.txt**  // The microbiome file which contains filtered, transformed and covariate-corrected microbial taxa. See [Chapter 1](#chapter-1-16s-data-processing) and [Step 4.1](#step-41-correct-microbiome-for-covariate-create-annotation-files)
2. **tax_numeric.txt.annot** //Annotation file for taxonomies. See [Step 4.1](#step-41-correct-microbiome-for-covariate-create-annotation-files)
3. **tax_binary.txt**  // Presence/absence microbiome file. See [Step 4.1](#step-41-correct-microbiome-for-covariate-create-annotation-files)
4. **tax_binary.txt.annot** //Presence/absence annotation file. See [Step 4.1](#step-41-correct-microbiome-for-covariate-create-annotation-files)
5. **trityper_chapter4** // The folder which contains MAF and QC filtered autosomal genotypes in TriTyper format. See [Step 4.2](#step-42-filter-genotypes)
6. **coupling_file.txt** // Coupling file thanks links microbiome IDs to genotype IDs. See [Step 4.3](#step-43-prepare-coupling-file)
5. **eqtl-mapping-pipeline-1.4nZ** // The folder with eQTL pipeline binaries
6. **settings.xml** // XML settings file for quantitative analysis. See [Step 4.4](#step-44-generate-xml-setting-files-for-quantitative-and-binary-analyses)
7. **settings_binary.xml** // XML settings file for binary analysis. See [Step 4.4](#step-44-generate-xml-setting-files-for-quantitative-and-binary-analyses)

If you have all of these ready, just run this commands from a *screen* session:

For quantitative trait analysis:
```
java -XX:ParallelGCThreads=5 -Xmx30G -jar eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings settings.xml
```

For binary trait analysis:
```
java -XX:ParallelGCThreads=5 -Xmx30G -jar eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings settings_binary.xml
```

Having this finished, you will generate two folders, **COHORT_full_quant** and **COHORT_full_bin** that should be transferred to UMCG servers, see [Chapter 7](#chapter-7-acquire-guest-account-and-upload-data). 

**Note: If you were involved in [Benchmark run](#chapter-3-perform-genome-wide-association-study-benchmark), you can use the same login name and ssh keys, it's not necessary to request access again.** 

# Chapter 5. run GWAS on alpha diversity

This chapter describes the process of running GWAS on microbiome alpha diveristy. The data required for this analysis was generated during previous steps. In particular:

1. **taxonomy_table.txt**. It was generated at [Step 1.2](#step-12-process-mapping-results)
2. **phenotypes.txt**, phenotype file with genetic principal components. The same as used at [Step 4.1](#step-41-correct-microbiome-for-covariates-create-annotation-files)
3. **genotypes** in TriTyper format. It was generated at the [Step 4.2](#step-42-filter-genotypes)
4. **linkage.txt**, phenotype-microbiome coupling file (optional). The same as was used at the [Step 4.1](#step-41-correct-microbiome-for-covariates-create-annotation-files)
5. **coupling_file.txt**, genotype-microbiome coupling file. It was  generated at [Step 4.3](#step-43-prepare-coupling-file)

## Step 5.1. Generate alpha diversity tables corrected for covariates 

In your project folder, run:

```
#If your sample ID's are matched between 16S and phenotype data
Rscript step5.1_create_alpha_diveristy.R taxonomy_table.txt phenotypes.txt

#If your sample ID's are different between phenotypes and 16s libraries
Rscript step5.1_create_alpha_diveristy.R taxonomy_table.txt phenotypes.txt linkage.txt
```

The files **alpha_div.txt** and **alpha_div.txt.annot** will be generated


## Step 5.2. Create XML settings file 

Replace MBRUN to the name of your cohort and run:

```
cat miQTL_cookbook/software/template_alpha.xml|perl -pe 's/COHORTNAME/MBRUN/' > ./settings_alpha.xml
```

## Step 5.3. Run alpha diversity GWAS

Please check that you have **trityper_chapter4** folder, **coupling_file.txt** file, and **eqtl-mapping-pipeline-1.4nZ** in the folder from where you are going to submit this task, see [Step 4.5](#step-45-run-full-miqtl-mapping) for details. After that, run:

```
java -XX:ParallelGCThreads=5 -Xmx30G -jar eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings settings_alpha.xml
```

The folder **COHORT_alpha** will be generated. Please upload it to our sftp server. See [Chapter 7](#chapter-7-acquire-guest-account-and-upload-data) for more details.

# Chapter 6. Run binary trait validation

To run microbiome-wide GWAS on bacterial presence, we used fast but not precised method. To obtain accurate significance estimation, it's necessary to re-analyze top hits with slow but more correct procedure. To do so, please do the following:

1. Update your local github repository copy and extract the list of SNP-bacterial pairs for further validaiton from [database folder](https://github.com/alexa-kur/miQTL_cookbook/tree/master/database):

```
zcat miQTL_Cookbook/database/topBinary_Sep2018.txt.gz > topBinary_Sep2018.txt
```

## Step 6.1 first batch of validation results

Check that you have the following files/folders in your project folder:
1. **tax_binary.txt** file. The one you generated at [Step 4.1](#step-41-correct-microbiome-for-covariates-create-annotation-files)
2. **phenotypes.txt** file. The one you used at [Step 4.1](#step-41-correct-microbiome-for-covariates-create-annotation-files)
3. **trityper_chapter4** genotypes folder. The one you generated at [Step 4.2](#step-42-filter-genotypes)
4. **GenotypeHarmonizer-1.4.20-SNAPSHOT** folder. You used it at [Step 4.2](#step-42-filter-genotypes)
5. **coupling_file.txt** file. The one you prepared at [Step 4.3](#step-43-prepare-coupling-file) 
6. **topBinary_Sep2018.txt**. See above. 

If necessary, you can use other names for these files/folders by using builtin script parameters: *-m, -p, -g, -h, -c and -r* respectively.


Then, run the script **step6.1_validationBinary.R**, replacing COHORT with a name of your cohort:

```
Rscript step6.1_validationBinary.R -n COHORT
```

It will generate the file named **COHORT_binaryResults_updated.txt.gz**. Please upload it to [SFTP](#chapter-7-acquire-guest-account-and-upload-data) or our [Dropbox](https://www.dropbox.com/home/Microbiome-QTL_Charge) in the folder **binary_validation** and send the notification to Alexander Kurilshikov (alexa.kur@gmail.com).


## Step 6.2 second batch of validation results

Since we got the updated summary statitistics after joing several new cohorts, it's necessary to re-evaluate several new SNP-bacterium associations. To do so, you need another script, **step6.2_validationBinary_step2.R**, and another database file **topBinary_Jan2019.txt**:

```
zcat miQTL_Cookbook/database/topBinary_Jan2019.txt.gz > topBinary_Jan2019.txt
```

The usage of the script is completely the same as previous one:

```
Rscript step6.2_validationBinary_step2.R -n COHORT
```

Again, you can can also modify the paths to the files (see [Step 6.1](#step-61-first-batch-of-validation-results)), if necessary. The output file called **COHORT_binaryResults_Jan2019_updated.txt.gz** should be sent to Alex Kurilshikov' e-mail (alexa.kur@gmail.com)




# Chapter 7. Acquire guest account and upload data 

1. Generate private and public key pair as documented here:
[http://wiki.gcc.rug.nl/wiki/RequestAccount]
2. Send the e-mail to GCC Helpdesk (Helpdesk dot GCC dot Groningen at gmail dot com) and Alex Kurilshikov (alexa.kur @ gmail.com) asking the guest account. The format of letter in the following:
```{r, eval = FALSE}
Dear GCC Helpdesk,

 I am participating MiBioGen Consortium,contact is Alex Kurilshikov (in cc). The name of the cohort is [LLD/MIBS/ROTTERDAM/FGFP, etc].

 Please provide me guest account to the cher-ami SFTP server for three months.
 
 My public key is attached.

 Please cc: alexa.kur@gmail.com and sasha.zhernakova@gmail.com to all the correspondence and provide them access to my guest account.

 Best regards,
 [Name]
 [Institution]
```

3. If you are ready with the analyses, please upload the the data to the server to your guest account acquired in Step 7 (named umcg-guest[0-9]). Please create the folder which is named by your cohort name and upload two output directories generated in Step 4, named **./COHORTNAME_full_quant** and **./COHORTNAME_full_bin**. You can use several options for uploading, like **lftp** or **sftp** 

```{r,eval = FALSE}
#start lftp

lftp
 
# connect with your guest account:

lftp :~> open -u [your_guest_accountname],none -p 22 sftp://cher-ami.hpc.rug.nl

# create remote directory named as your cohort, say, in this case it will be Rotterdam 

lftp :~> mkdir Rotterdam

# go to this directory

lftp :~> cd Rotterdam

# copy your output for quantitative trait analysis, assuming it's named as Rotterdam_output

lftp :~> mirror -R Rotterdam_output

# copy your output for binary trait analysis, assuming it's named as Rotterdam_output_binary

lftp :~> mirror -R Rotterdam_output_binary

# exit from remote server

lftp :~> exit
```


Now it's our job to run meta-analysis.




# Additions

## Appendix 1. OTU-based 16S data processing

For those who are still interested in performing 16S processing for their own purposes, we provide the OTU-based pipeline. Please note that for miQTL meta-analysis we require using RDPclassifer pipeline described in section 1!

### Installation

First we need to install QIIME on your machine. I should be any kind of Linux machine, including physical or virtual machine. You don't need to have root priviliges to use QIIME

```{bash,eval = FALSE}
cd ~/
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod a+x Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
#follow installation instructions and install Miniconda in ~/miniconda3 folder
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
source activate qiime1
#check if installation was finished succesfully
print_qiime_config.py
```

### SILVA database downloading

```{bash,eval = FALSE}
#create project dir
mkdir 16S_picking
#copy your 16S data in the folder
cp PATH_TO_YOUR_FASTA_FILE ./16S_picking
cd 16S_picking

#download SILVA v.119 database
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_119_release.zip
unzip Silva_119_release.zip
#create custom config file
echo "\
pick_otus_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
pynast_template_alignment_fp $PWD/Silva119_release/core_alignment/core_Silva119_alignment.fna
assign_taxonomy_reference_seqs_fp $PWD/Silva119_release/rep_set/97/Silva_119_rep_set97.fna
assign_taxonomy_id_to_taxonomy_fp $PWD/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt" > qiime_config
```

### Running the OTU picking (closed OTU picking method with 97% cutoff)

```{bash,eval = FALSE}
# if you have multiple cores, you can run analysis faster in parallel, redefining the $PROCNUM varialbe 
PROCNUM=1
source activate qiime1
export QIIME_CONFIG_FP=$PWD/qiime_config
pick_closed_reference_otus.py -i FASTA_FILE -o RESULT_FOLDER -a -O $PROCNUM
cd RESULT_FOLDER
biom convert -i otu_table.biom -o otu_table.tsv --to-tsv --header-key taxonomy
cat otu_table.tsv|tail -n+2 |perl -pe "s/#OTU ID/OTU_ID/" > temp.tsv
mv temp.tsv otu_table.tsv
```

### Getting taxonomies from OTU table

For this step, you should have R to be installed. 

```{r, eval = FALSE}
cutoff_presence = 0.1
get_taxonomy_table = function(otu_table, replace_string,cutoff = 0.1){
  otu_notax = as.matrix(otu_table[,-ncol(otu_table)])
  taxonomy = sub(replace_string,"",otu_table[,ncol(otu_table)])
  dnew = aggregate(otu_notax ~ as.factor(taxonomy),FUN = sum)
  rownames(dnew) = as.character(dnew[,1])
  dnew = dnew[,-1]
  dnew = t(dnew)
  dnew
  filter = dnew[,(colSums(dnew > 0) > 0.1*nrow(dnew))]
  return(dnew)
}
otus = read.table("otu_table.tsv",header=T,row.names=1,sep="\t",as.is = T)
metadata = data.frame(tax = c("genus","family","order","class"),replace_string = c("; D_6.*","; D_5.*","; D_4.*","; D_3.*"))
result = list()
for (i in 1:nrow(metadata)){
  taxonomy_table = get_taxonomy_table(otus,replace_string = metadata[i,2])
  result[[i]] = taxonomy_table
}
final_table = cbind(result[[1]],result[[2]],result[[3]],result[[4]])
final_table = final_table[,(colSums(final_table>0)>cutoff_presence * nrow(final_table))]
final_ln = log(final_table)
final_table[final_table == "-Inf"] = NA
write.table(final_table,file = "microbes.txt",sep="\t")
```

## Appendix 2. FASTA file rarefaction

We recommend to perform rarefaction before running the pipeline. To make it truly random, this procedure can be applied:

1. Install QIIME (see Appendix 1)
2. Create merged file for all samples with number of reads larger than threshold you want 
2. Check if your fasta headers are formatted properly (see section 1)
2. Download script **run_rarefaction.R** from Consofrtium dropbox [CookBook_CHARGE_miQTL/software](https://www.dropbox.com/home/Microbiome-QTL_Charge/CookBook_CHARGE_miQTL/software)
3. Put it into your project folder. Replace SEQUENCES.FASTA to the file where you store your sequences and run:

```{}
Rscript step0.2_run_rarefaction.R SEQUENCES.FASTA 10000
source activate qiime1
filter_fasta.py -f SEQUENCES.FASTA -o SEQ_RARIFIED.FASTA -s 2filter.ids
rm 2filter.ids

```

It will generate rarefied fasta file called **SEQ_RARIFIED.FASTA**

