
rm(list = ls())

library(broom)

## -----------------------------------------------

## global settings
if(TRUE) {
  
  cohorts = c('ATM_PAAD','TCGA_PAAD')
  # cohort = 'TCGA_PAAD' ## ATM_PAAD TCGA_PAAD
  altDP = 5
  
}

## -----------------------------------------------

for(my.analysis in 1:3) {
  print(my.analysis)
  funcs = c()
  if(my.analysis == 1) {
    funcs = c('De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins','Start_Codon_SNP')
  } else if (my.analysis ==2) {
    funcs = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins')
  } else if (my.analysis ==3) {
    funcs = c('Frame_Shift_Del','Frame_Shift_Ins','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins')
  }

  ## -----------------------------------------------
  
  ## cohorts 
  if (TRUE) { 
    for(cohort in cohorts) {
      
      sm.total = 0
      if(cohort == 'ATM_PAAD') { sm.total = 21}
      if(cohort == 'TCGA_PAAD') { sm.total = 150}
      
      
      ## flags
      ## whether to include cnv data into summary stats 
      for(use.cnv in c(0,1)) {
        
        cnv.outstring = ''
        if(use.cnv==1) { cnv.outstring='.wCNV'}
        
        for(ex.TP53 in c(0, 1)) { 
          
          TP53.outstring = ''
          if(ex.TP53 == 1) { TP53.outstring = '.exTP53'}
          
          ## -----------------------------------------------
            
          rm(maf.file, maf, cnv.file, cnv)
            
          ## -----------------------------------------------
          
          path = 'results/exome_genome'
          setwd(path)
          
          ## -----------------------------------------------
          
          ## import maf files 
          if(TRUE) {
          
            if(cohort == 'ATM_PAAD') {
              
              ## couch sm21
              maf.dir = path
              maf.file = paste0('ATM_PAAD.exome_genome.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP',altDP,'.NSSM.add4Genes.maf')
              
              
            } else if (cohort == 'TCGA_PAAD') {
              ## tcga prad
              maf.dir = '/01_Data/GDAC_data/mutsig'
              maf.file = 'PAAD-TP.final_analysis_set.maf.fixed.slim.maf'
          
              
            }
            
            print(maf.file)
            
            ## -----------------------------------------------
            
            maf = read.delim(file.path(maf.dir,maf.file),stringsAsFactors = F)
            dim(maf) ## 2246  282
            
            print(dim(maf[maf$Hugo_Symbol=='TP53',,drop=F]))
            
            table(maf$Variant_Classification)
            
            dim(maf) ## 2246  282
            maf = maf[!is.na(maf$Variant_Classification) &
                        maf$Variant_Classification %in% funcs,,drop=F]
            dim(maf) ## 2238  282
            
            ## -----------------------------------------------
            
            if(cohort == 'TCGA_PAAD') {
              
              tcga.paad.clinical = read.delim('sampleinfo/TCGA_PAAD.Table_S1.txt',stringsAsFactors = F)
      
              dim(tcga.paad.clinical) ## 150  53
              
              ## -----------------------------------------------
              
              maf  = maf[gsub('[-]\\w{3}[-]\\w{4}[-]\\w{2}$','',
                              maf$Tumor_Sample_Barcode) %in%tcga.paad.clinical$Tumor.Sample.ID,]
              dim(maf) ## 40505    37
              
              atm.samples = tcga.paad.clinical$Tumor.Sample.ID[
                grep('ATM',tcga.paad.clinical$Mutated.PanCan.Genes)
                ]
              
              print(atm.samples)
            }
            
          
          }
          
          ## -----------------------------------------------
          
          if(use.cnv == 1) {
            
            if(cohort == 'ATM_PAAD') {
              
              ## couch sm21
              cnv.file = 'results/genome/gistic2/all_thresholded.by_genes.txt'
              
              cnv = read.delim(cnv.file, stringsAsFactors = F)
              
            } else if (cohort == 'TCGA_PAAD') {
              
              ## tcga sm150
              load('/01_Data/GDAC_data/cnv/TCGA_allCancer.copynumber.all_thresholded_by_gene.txt.RData')
              dim(data)
              cnv = data
              rm(data)
              
              cnv = cnv[,gsub('[.]\\w{3}[.]\\w{4}[.]\\w{2}$','',colnames(cnv)) %in%
                          c('Gene.Symbol','Locus.ID','Cytoband',
                            gsub('[.]\\w{3}[.]\\w{4}[.]\\w{2}$','',
                                 gsub('-','.',maf$Tumor_Sample_Barcode)))]
          
            }
            
            ## -----------------------------------------------
            
            print(dim(cnv)) ## 25988    18
            row.names(cnv) = cnv$Gene.Symbol
            cnv = cnv[,-c(1:3)]
            
            print(dim(cnv)) ## 25988    15
            
          }
          
          ## -----------------------------------------------
          
          if(ex.TP53==1) {
            maf = maf[!maf$Hugo_Symbol %in% c('TP53'),]
            if(use.cnv == 1) { cnv = cnv[!row.names(cnv) %in% c('TP53'),] }
          }
          
          ## -----------------------------------------------
          
          ## compute mut prevelance 
          if(TRUE) {
            
            ## mutated gene list
            if(TRUE) {
              my.df = NULL 
              my.df = maf[,c('Hugo_Symbol','Tumor_Sample_Barcode')]
              print(dim(my.df))
              
              if(use.cnv==1) {
                my.df.2 = NULL 
                my.df.2 = cnv
                
                if(TRUE) {
                  my.df.3 = NULL 
                  my.df.3 = reshape2::melt(data.frame(
                    Gene = row.names(my.df.2),
                    my.df.2,
                    stringsAsFactors = F
                  ))
                  
                  ## only keep deep deletions!!
                  my.df.3 = my.df.3[my.df.3$value == -2,,drop=F]
                  colnames(my.df.3) = c('Hugo_Symbol','Tumor_Sample_Barcode','cn')
                  
                  my.df.3$Tumor_Sample_Barcode = gsub('[.]','-',my.df.3$Tumor_Sample_Barcode) 
                }
                
                my.df = rbind(
                  my.df[,c('Hugo_Symbol','Tumor_Sample_Barcode')],
                  my.df.3[,c('Hugo_Symbol','Tumor_Sample_Barcode')]
                  
                )  
              }
              
              if(cohort =='TCGA_PAAD') {
                my.df$Tumor_Sample_Barcode = gsub('[-]\\w{3}[-]\\w{4}[-]\\w{2}$','',
                                                  my.df$Tumor_Sample_Barcode)
              }
          
              my.df = unique(my.df)
              my.df = my.df[order(my.df$Hugo_Symbol),]
              
              my.df = data.frame(table(my.df$Hugo_Symbol))
              colnames(my.df) = c('Gene','samples.total')
              
              my.df$samples.sequenced =sm.total
              my.df$samples.frac =  my.df$samples.total / my.df$samples.sequenced 
            }
            
            ## -----------------------------------------------
            
            write.table(my.df,
                      file = paste0(cohort,
                                    cnv.outstring,TP53.outstring,
                                    '.func',length(funcs),'.',
                                    maf.file,'.genelist.txt'),
                      sep = '\t', quote = F,row.names = F)
            
          }
            
        }
      
      }
      
    }
  }
  
}

## -----------------------------------------------

sessionInfo()

