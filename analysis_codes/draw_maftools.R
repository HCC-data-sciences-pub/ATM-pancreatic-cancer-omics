
rm(list = ls())

if(TRUE) {
  library(cluster)
  library(Biobase)
  library(ctc)
  library(ape)
  library(ggdendro)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(gplots)
  library(stringr)
  library(RColorBrewer)
  
  library(ComplexHeatmap)
  library(broom)
  library(survival)
  library(survminer)
  library(circlize)
  
  ## https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
  ## https://github.com/cran/ggsci/blob/master/R/discrete-nejm.R
  library(ggsci)
  library(scales)
  library(ggpubr)
  library(plotrix)
  library(maftools)
  library(edgeR)
  library(xCell)
  library(fBasics)
}

## --------------------------------------------------------

## settings
if(TRUE) {
  
  rnaseq.sample.ex = c('T7290343', '5436736_Repeat')
  exome.sample.ex = c()
  sample.order = c('PNATM-123-tumor', 'PNATM-108-tumor', 'PNATM-112-tumor', 'PNATM-119-tumor', 'PNATM-113-tumor', 'PNATM-109-tumor', 'PNATM-101-tumor', 'PNATM-102-tumor', 'PNATM-104-tumor', 'PNATM-105-tumor', 'PNATM-106-tumor', 'PNATM-114-tumor', 'PNATM-115-tumor', 'PNATM-116-tumor', 'PNATM-117-tumor', 'PNATM-118-tumor', 'PNATM-121-tumor', 'PNATM-122-tumor', 'PNATM-125-tumor', 'PNATM-103-tumor', 'PNATM-124-tumor')
  
  output = 'ATM_PAAD.exome_genome'
  
}

## --------------------------------------------------------

## data files 
if(TRUE) {
  
  clinical.file = '../../sampleinfo/Sample list 2.0.RB.txt'
  maf.file = paste0(output,'.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP5.NSSM.add4Genes.maf')

  
}

## --------------------------------------------------------

path = 'results/exome_genome'
setwd(path)

## --------------------------------------------------------

## only run once: filter variants 
if(FALSE) {
  
  clinical = read.delim('../../sampleinfo/Sample list 2.0.RB.txt', stringsAsFactors = F)
  dim(clinical) ## 25
  mafs = NULL 
  
  ## import 
  for(my.sm in sort(clinical$SampleID.DNAseq)) {
    
    print(my.sm)
    my.df = NULL 
    my.file = NULL 
    
    my.file = paste0('funcotator/',my.sm,'-tumor.bwamem.mutect2.flt.funcotator.fixed.maf')
    if(! file.exists(my.file)) { next }
    
    my.df = read.delim(my.file, stringsAsFactors = F, skip = 2852)
    print(dim(my.df))
    
    mafs = rbind(mafs,
                 my.df)
    
  }
  
  print(dim(mafs)) ## 13894   301
  table(mafs$Tumor_Sample_Barcode)
  
  write.table(mafs,
              file = paste0(output,'.maf'), sep = '\t', row.names = F, quote = F)
  
  ## --------------------------------------------------------
  
  ## filter 
  dim(mafs) ## 14121   301

  mafs = read.delim(paste0(output,'.fixed.maf'), stringsAsFactors = F,
                    na.strings = c('NA','na','','N/A','-','.'))
  class(mafs$gnomAD_exome_AF)

  mafs = mafs[is.na(mafs$gnomAD_exome_AF) | 
                (!is.na(mafs$gnomAD_exome_AF) & (
                  mafs$gnomAD_exome_AF == '' | 
                    mafs$gnomAD_exome_AF < 0.0001
                )),]
  dim(mafs) ## 12588   301
  
  write.table(mafs,
              file = paste0(output,'.fixed.gnomAD_exome_AF0.0001.maf'), sep = '\t', row.names = F, quote = F)

  mafs.0 = mafs
  basicStats(mafs.0$DP)
  # X..mafs.0.DP
  # nobs        1.258800e+04
  # NAs         0.000000e+00
  # Minimum     2.300000e+01
  # Maximum     3.293000e+03
  # 1. Quartile 1.147500e+02
  # 3. Quartile 2.770000e+02
  # Mean        2.201122e+02
  # Median      1.780000e+02
  # Sum         2.770773e+06
  # SE Mean     1.448447e+00
  # LCL Mean    2.172731e+02
  # UCL Mean    2.229514e+02
  # Variance    2.640961e+04
  # Stdev       1.625103e+02
  # Skewness    3.587103e+00
  # Kurtosis    3.120800e+01
  
  mafs = mafs.0[ !is.na(mafs.0$t_alt_count) &
                   mafs.0$t_alt_count>=8,]
  
  dim(mafs) ## 2621  301
  
  write.table(mafs,
              file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altDP5.maf'), sep = '\t', row.names = F, quote = F)
  
  mafs = mafs.0[!is.na(mafs.0$tumor_f) &
                  mafs.0$tumor_f >= 0.1 & 
                !is.na(mafs.0$t_alt_count) &
                  mafs.0$t_alt_count>=4,]
  
  dim(mafs) ## 1609  301
  
  write.table(mafs,
              file = paste0(output,'.fixed.gnomAD_exome_AF0.0001_altAF0.1_altDP4.maf'), sep = '\t', row.names = F, quote = F)
  
}

## --------------------------------------------------------

## import data 
if(TRUE) {
  
  clinical = read.delim(clinical.file, stringsAsFactors = F)
  clinical$Tumor_Sample_Barcode = paste0(clinical$SampleID.DNAseq,'-tumor')
  
  maf = read.maf(maf.file, clinicalData = clinical)
  
}

## --------------------------------------------------------

## prase maf
if(TRUE) {
  
  my.genes = c('TP53','KRAS','CDKN2A','SMAD4')
  
  maf.sampleSummary = getSampleSummary(maf)
  
  clinical[!clinical$Tumor_Sample_Barcode %in% maf.sampleSummary$Tumor_Sample_Barcode,]
  
  maf.geneSummary = getGeneSummary(maf)
  maf.geneSummary[maf.geneSummary$Hugo_Symbol %in% my.genes,]
  
  ## --------------------------------------------------------
  
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', 
                 dashboard = TRUE, titvRaw = FALSE)
  
  ## oncoplot for top ten mutated genes.
  oncoplot(maf = maf, top = 30, removeNonMutated = F,
           clinicalFeatures = c('Stage.at.diagnosis'))
  
  ## --------------------------------------------------------
  
  ## 2 oncoprints
  pdf(file = paste0(output,'.oncoprint.4genes.pdf'),
      width = 6, height = 4)
  oncoplot(maf = maf,
           genes = my.genes,
           sampleOrder = sample.order,
           removeNonMutated = F,
           showTumorSampleBarcodes = T,
           clinicalFeatures = c('Gender','Stage.at.diagnosis'))
  dev.off()
  ## manually retrieve sample order here ....
  
  pdf(file = paste0(output,'.oncoprint.pdf'),
      width = 10, height = 20)
  oncoplot(maf = maf,
           top = 60,
           sampleOrder = sample.order,
           removeNonMutated = F,
           showTumorSampleBarcodes = T,
           clinicalFeatures = c('Stage.at.diagnosis'))
  dev.off()
  
  
  
  maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
  #plot titv summary
  plotTiTv(res = maf.titv)
  
  plotVaf(maf = maf, vafCol = 'tumor_f')
  
  #exclusive/co-occurance event analysis on top 10 mutated genes. 
  somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
  
  ## Detecting cancer driver genes based on positional clustering
  maf.sig = oncodrive(maf = maf, AACol = 'Protein_Change', 
                      minMut = 5, pvalMethod = 'zscore')
  head(maf.sig)
  
  plotOncodrive(res = maf.sig, fdrCutOff = 0.15, useFraction = TRUE)
  rm(maf.sig)
  
  lollipopPlot(maf = maf, 
               gene = 'MTCH2', labelPos = 'all',
               AACol = 'Protein_Change', showMutationRate = TRUE)
  
  lollipopPlot(maf = maf, 
               gene = 'KRAS', labelPos = 'all',
               AACol = 'Protein_Change', showMutationRate = TRUE)
  
  ## Adding and summarizing pfam domains
  pfamDomains(maf = maf, AACol = 'Protein_Change', top = 10)
  
  ## total TMB 
  maf.samplesummary = getSampleSummary(maf)
  clinical$Tumor_Sample_Barcode = as.character(clinical$Tumor_Sample_Barcode)
  maf.samplesummary$Tumor_Sample_Barcode = as.character(maf.samplesummary$Tumor_Sample_Barcode)
  
  clinical = merge(clinical, maf.samplesummary, by = 'Tumor_Sample_Barcode')
  
  write.csv(clinical,
            file = paste0(clinical.file,'.wSeqfiles.wTMB.csv'))
  
  p1 = ggplot(clinical, aes(Stage.at.diagnosis, total)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
    # scale_fill_manual(values = plot.colors) + 
    # stat_compare_means(method = 'wilcox.test') +
    ggtitle('Total TMB (Non-Synonymous Somatic Mutations (SNVs/indels), NSSMs)') +
    theme_classic()
  
  p1
  
  p2 = ggplot(clinical, aes(Location, total)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
    # scale_fill_manual(values = plot.colors) + 
    stat_compare_means(method = 'wilcox.test') +
    ggtitle('Total TMB (Non-Synonymous Somatic Mutations (SNVs/indels), NSSMs)') +
    theme_pubr()
  
  p2
  
  ## log transform ...
  min(clinical$total) ## 2
  max(clinical$total) ## 337
  dim(clinical) ## 40 samples 
  
  p2.1 = ggplot(clinical, aes(Group.2, total)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
    geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
    # scale_fill_manual(values = plot.colors) + 
    stat_compare_means(method = 't.test') +
    ggtitle('Total TMB (Non-Synonymous Somatic Mutations (SNVs/indels), NSSMs)') +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'l')  +
    theme_pubr()
  
  p2.1
  
  p3 = ggplot(clinical, aes(total)) +
    geom_density()
  
  p3
  
  p3.1 = ggplot(clinical, aes(log10(total))) +
    geom_density()
  
  p3.1
  
  
}

## --------------------------------------------------------

## plot KRAS
if(TRUE) {
  
  maf.file.0 = paste0(output,'.fixed.maf')
  purity.file = 'ATM_PAAD.bwamem.merged.rmdup.bam_info.txt'
  kras.steven.file = 'kras_mutations.steven.txt'
  
  maf.0 = read.maf(maf.file.0, clinicalData = clinical)
  purity = read.delim(purity.file, stringsAsFactors = F, header = F)
  kras.steven = read.delim(kras.steven.file, stringsAsFactors = F)
  
  ## --------------------------------------------------------
  
  maf.sampleSummary = getSampleSummary(maf.0)
  clinical[!clinical$Tumor_Sample_Barcode %in% maf.sampleSummary$Tumor_Sample_Barcode,]
  
  
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', 
                 dashboard = TRUE, titvRaw = FALSE)
  
  ## oncoplot for top ten mutated genes.
  oncoplot(maf = maf.0, top = 30, removeNonMutated = F,
           clinicalFeatures = c('Stage.at.diagnosis'))
  
  maf.0.sub = subsetMaf(maf.0, genes = c('KRAS'))
  
  getSampleSummary(maf.0.sub)
  
  write.mafSummary(maf.0.sub, basename = paste0(output,'.fixed.kras'))
  
  ## --------------------------------------------------------
  
  purity$Sample = gsub('-tumor.+$','',purity$V1)
  
  purity = reshape2::dcast(Sample ~ V2, data = purity, value.var = 'V3')
  
  write.csv(purity,
            file = paste0(purity.file,'.matrix.csv'))
  
  ## --------------------------------------------------------
  
  sommut = read.delim('ATM_PAAD.exome.fixed.kras_maftools.maf', stringsAsFactors = F)
  
  sommut$Sample = gsub('-tumor','',sommut$Tumor_Sample_Barcode)
  
  sommut = sommut[,c('Sample','Hugo_Symbol','Protein_Change','tumor_f', 't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count')]
  
  # sommut = merge(sommut, purity[,c('Sample','Output_Ploidy','Sample_Purity')])
  
  kras.steven$Sample = gsub('^PTATM-0','PNATM-1',kras.steven$Sample)
  
  kras.steven = merge(kras.steven, sommut, by = 'Sample', all = T)
  
  kras.steven = merge(kras.steven, purity[,c('Sample','Output_Ploidy','Sample_Purity')], all = T)
  
  write.csv(kras.steven, 
            file = paste0(kras.steven.file,'.wSommut.csv'))
  
  ## --------------------------------------------------------
  
  data.plot = kras.steven
  data.plot$group = NA
  data.plot$group[which(data.plot$X..Reads.Supporting>=20)] = 'high'
  data.plot$group[which(data.plot$X..Reads.Supporting<20)] = 'low'
  
  
  data.plot$Sample_Purity = as.numeric(as.character(data.plot$Sample_Purity))
  
  p1 = ggplot(data.plot, aes(group, Sample_Purity)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 2) +
    geom_jitter(width = 0.1, height = 0, size = 2) +
    theme_pubr() +
    stat_compare_means(method = 'wilcox.test')

  p1  
  
  
  p1 = ggplot(data.plot, aes(group, t_alt_count)) +
    geom_boxplot(width = 0.5, color = '#C0C0C0', size = 2) +
    geom_jitter(width = 0.1, height = 0, size = 2) +
    theme_pubr() +
    stat_compare_means(method = 'wilcox.test')
  
  p1  
  
  
}


## --------------------------------------------------------


sessionInfo()
