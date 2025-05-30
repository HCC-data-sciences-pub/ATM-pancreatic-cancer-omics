

rm(list = ls())

if(TRUE) {
  
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(broom)
  library(fBasics)
  library(UpSetR)
  library(ComplexHeatmap)
  # library(VariantAnnotation)
  library(sigminer)
  # library(NMF)
  library(maftools)
  
}

## --------------------------------------------------

path = 'results/exome_genome'
setwd(path)

## --------------------------------------------------

## only run once 
## merge genome maf and exome maf
if(FALSE) {

  maf.exome.file = '../exome/ATM_PAAD.exome.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP8.maf'
  maf.genome.file = '../genome/ATM_PAAD.genome.bwamem.mutect2.gnomAD_genome_AF0.0001_altDP8.maf'
  
  maf.exome = read.delim(maf.exome.file, stringsAsFactors = F)
  maf.genome = read.delim(maf.genome.file, stringsAsFactors = F)
  
  table(maf.exome$Source_MAF)
  table(maf.genome$Source_MAF)
  
  maf.exome$Key = paste0(maf.exome$Tumor_Sample_Barcode,'!',
                         maf.exome$Chromosome,'!',
                         maf.exome$Start_Position,'!',
                         maf.exome$End_Position,'!',
                         maf.exome$Reference_Allele,'!',
                         maf.exome$Tumor_Seq_Allele2)
  maf.genome$Key = paste0(maf.genome$Tumor_Sample_Barcode,'!',
                          maf.genome$Chromosome,'!',
                          maf.genome$Start_Position,'!',
                          maf.genome$End_Position,'!',
                          maf.genome$Reference_Allele,'!',
                          maf.genome$Tumor_Seq_Allele2)
  
  mutations = sort(unique(c(maf.exome$Key, maf.genome$Key)))
  
  length(mutations) ## 30828
  
  maf = NULL 
  for(i in 1:length(mutations)) {
    # print(i)
    
    if (i %% 1000 == 0) { print(i)}
    
    my.mut = mutations[i]
    
    my.df = NULL 
    
    if(my.mut %in% maf.exome$Key) {
      my.df = maf.exome[maf.exome$Key == my.mut,,drop=F]
    } else if (my.mut %in% maf.genome$Key) {
      my.df = maf.genome[maf.genome$Key == my.mut,,drop=F]
    } else {
      print(paste0("mutation does not exist!", my.mut))
    }
    
    maf = rbind(maf,
                my.df)
    
  }
  
  maf = maf[!colnames(maf) %in% c('Key'),]
  
  write.table(maf,
              file = 'ATM_PAAD.exome_genome.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP8.maf',
              sep = '\t', row.names = F, quote = F)
  
  
}

## --------------------------------------------------

## import data
if(TRUE) {
  
  # maf = read.delim(maf.file, stringsAsFactors = F)
  # dim(maf) ## 2621  307
  
  maf =  read_maf(maf = maf.file)
}