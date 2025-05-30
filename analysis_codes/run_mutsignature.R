

rm(list = ls())

if(TRUE) {
  
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(broom)
  library(fBasics)
  library(UpSetR)
  library(ComplexHeatmap)
  library(VariantAnnotation)
  library(sigminer)
  
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  
}

## --------------------------------------------------

path = 'results'
setwd(path)

## --------------------------------------------------

maf.file = 'exome/ATM_PAAD.exome.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP8.maf'
maf.file = 'genome/ATM_PAAD.genome.bwamem.mutect2.gnomAD_genome_AF0.0001_altDP8.maf'
vcf.file = paste0(maf.file, '.vcf')

## --------------------------------------------------

## only run once 
## split maf into per sample maf!
if(FALSE) {
  
  maf = read.delim(maf.file, stringsAsFactors = F)
  
  samples = sort(unique(gsub('.bwamem.mutect2.flt.funcotator.fixed','',maf$Source_MAF)))
  
  for(my.sm in samples) {
    
    print(my.sm)
    
    my.df = NULL 
    my.df = maf[gsub('.bwamem.mutect2.flt.funcotator.fixed','',maf$Source_MAF) == my.sm,,drop=F]
    
    write.table(my.df,
                file = paste0(maf.file,'.',my.sm,'.maf'),
                sep = '\t', col.names = T, row.names = F, quote = F)
    
  }
  
  
}

## --------------------------------------------------

for(my.sm in samples) {
  
  print(my.sm)
  maf.file = paste0(maf.file,'.',my.sm,'.maf')
  
  

## --------------------------------------------------

## import data
if(TRUE) {
  
  # maf = read.delim(maf.file, stringsAsFactors = F)
  # dim(maf) ## 2621  307
  
  maf =  read_maf(maf = maf.file)
}

## --------------------------------------------------

## de novo detection of signature
if(TRUE) {
  
  head(maf@data)
  
  head(maf@maf.silent)
  
  slotNames(maf)
  
  mt_tally <- sig_tally(
    maf,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
    useSyn = TRUE
  )
  
  mt_tally$nmf_matrix[1:5, 1:5]
  
  samples=sort(unique(row.names( mt_tally$nmf_matrix)))
  
  for(my.sm in samples) {
    print(my.sm)
    
    p1 = NULL 
    p1 = show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", x_label_angle = 90,
                   samples_name = my.sm,
                   samples = my.sm)
    pdf(file = paste0(maf.file,'.',my.sm,'.sigminer.raw.pdf'),
        height = 3, width = 10)
    print(p1)
    dev.off()
  }
  
  
  
  # ## --------------------------------------------------
  # 
  # mt_tally_all <- sig_tally(
  #   maf,
  #   ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  #   useSyn = TRUE,
  #   mode = "ALL",
  #   add_trans_bias = TRUE
  # )
  # 
  # str(mt_tally_all, max.level = 1)
  
}

## --------------------------------------------------

## extract signature 
if(FALSE) {
  
  mt_est <- sig_estimate(mt_tally$nmf_matrix,
                         range = 2:10,
                         nrun = 10, # increase this value if you wana a more stable estimation
                         use_random = FALSE, # if TRUE, add results from randomized input
                         cores = 4,
                         pConstant = 1e-13,
                         verbose = TRUE
  )


  ## You can also select the measures to show
  ## by 'what' option
  show_sig_number_survey2(mt_est$survey)

  show_sig_number_survey(mt_est$survey)

  ## --------------------------------------------------
  
  mt_sig <- sig_extract(mt_tally$nmf_matrix,
                        n_sig = 6,
                        nrun = 30,
                        cores = 4,
                        pConstant = 1e-13
  )
  
  sim_v3 <- get_sig_similarity(mt_sig, sig_db = "SBS")
  
  pheatmap::pheatmap(sim_v3$similarity)
  
  ## --------------------------------------------------
  
  sig_signature(mt_sig)[1:5, ]
  
  sig_exposure(mt_sig)[, 1:5]
  
  get_sig_exposure(mt_sig)
  
}

## --------------------------------------------------

## make plots
if(FALSE) {
  
  show_catalogue(t(mt_tally$nmf_matrix), style = "cosmic", x_label_angle = 90,
                 samples_name = my.sm,
                 samples = my.sm)
  
  ## --------------------------------------------------
  
  show_sig_profile(mt_sig, mode = "SBS", style = "cosmic", x_label_angle = 90)
  
  show_cosmic_sig_profile(sig_index = c(1, 6, 87,15), style = "cosmic")
  
  ## --------------------------------------------------
  
  show_sig_consensusmap(mt_sig)
  
}

## --------------------------------------------------

## convert maf back to vcf; only run once 
if(FALSE) {

  vcf = data.frame(CHROM = maf$Chromosome,
                   POS = maf$Start_Position,
                   ID = '.',
                   REF = maf$Reference_Allele,
                   ALT = maf$Tumor_Seq_Allele2,
                   QUAL = 60,
                   FILTER = 'PASS',
                   INFO = '.',
                   FORMAT = '.',
                   stringsAsFactors = F)
  
  dim(vcf) ## 2621    9
  
  colnames(vcf)[1] = '#CHROM'
  
  write.table(vcf,
              file = paste0(maf.file,'.vcf.body'),
              sep = '\t', row.names = F, quote = F)
  
}

## --------------------------------------------------

## mutation signatures 2-letter
if(FALSE) {
  
  data.plot = maf
  data.plot = data.plot[intersect(grep('^[ATGC]$',maf$Reference_Allele) , 
                                  grep('^[ATGC]$',data.plot$Tumor_Seq_Allele2)),]
  dim(data.plot) ## 2203  307
  
  data.plot$key = paste0(data.plot$Reference_Allele,':',
                         data.plot$Tumor_Seq_Allele2)
  
  table(data.plot$key)
  
  data.plot$Sample = data.plot$Tumor_Sample_Barcode
  
  data.plot = data.frame(table(data.plot[,c('Sample','key')]))
  data.plot
  
  # data.plot$Sample = factor(data.plot$Sample, levels = plot.order)
  
  p1 = ggplot(data.plot, aes(key, Freq)) +
    geom_bar(aes(fill = Sample), stat = 'identity') +
    # scale_fill_jco() +
    ggtitle('variants PASS, before variant function filter') +
    theme_pubr() +
    facet_grid(~Sample) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  print(p1)
  
  write.csv(data.plot,
            file = paste0(maf.file,
                          '.mutsig_2letter.bar.csv'))
  
  pdf(file = paste0(maf.file,
                    '.mutsig_2letter.bar.pdf'),
      width = 10, height = 5)
  print(p1)
  dev.off()
  
}


}

sessionInfo()

