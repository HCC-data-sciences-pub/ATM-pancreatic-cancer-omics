
rm(list = ls())

## reduce altDP filter from 8 to 5
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
  library(data.table)
  library(ggrepel)
  library(ggfortify)
}

## --------------------------------------------------------

## settings
if(TRUE) {
  
  # rnaseq.sample.ex = c('T7290343', 'T5436736_Repeat')
  ## still include T5436736_Repeat .. becasue it is ARID1A mutant!!!!!
  rnaseq.sample.ex = c('T7290343')
  exome.sample.ex = c()
  
  cohort = 'ATM_PAAD'
  # output = 'ATM_PAAD.exome_genome' ## ATM_PAAD.exome ATM_PAAD.exome
  altDP = 5
  
  my.analysis = 1
  use.cnv = 0
  ex.TP53 = 0 # default is 1; changed 1 to 0 for TP53 mutant vs nonmutant expr analysis 
  
  funcs = c()
  if(my.analysis == 1) {
    funcs = c('De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins','Start_Codon_SNP')
  } else if (my.analysis ==2) {
    funcs = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins')
  } else if (my.analysis ==3) {
    funcs = c('Frame_Shift_Del','Frame_Shift_Ins','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Start_Codon_Del','Start_Codon_Ins')
  }
  
  cnv.outstring = ''
  if(use.cnv==1) { cnv.outstring='.wCNV'}
  
  TP53.outstring = ''
  if(ex.TP53 == 1) { TP53.outstring = '.exTP53'}
  
  output = paste0(cohort,
                  cnv.outstring,TP53.outstring,
                  '.func',length(funcs))
  
  
  plot.colors = c('ATM_PAAD' = '#000000','TCGA_PAAD' = '#C0C0C0',
                 'mut' = '#CC0000','nonmut' = '#00CC00')
}

## --------------------------------------------------------

path = 'results/exome_genome'
setwd(path)

## --------------------------------------------------------

## data files and import data 
if(TRUE) {
  
  clinical.file = '../../sampleinfo/Sample list 2.0.RB.txt'
  clinical = read.delim(clinical.file, stringsAsFactors = F)
  clinical$Tumor_Sample_Barcode = paste0(clinical$SampleID.DNAseq,'-tumor')
  
  ## -----------------------------------------------
  
  ## import maf files 
  if(TRUE) {
    
    # if(cohort == 'ATM_PAAD') {
      
      ## couch sm21
      maf.dir = path
      maf.file = paste0('ATM_PAAD.exome_genome.bwamem.mutect2.gnomAD_exome_AF0.0001_altDP',altDP,'.NSSM.add4Genes.maf')
      
      
    # } else if (cohort == 'TCGA_PAAD') {
      ## tcga prad
      tcga.maf.dir = '01_Data/GDAC_data/mutsig'
      tcga.maf.file = 'PAAD-TP.final_analysis_set.maf.fixed.slim.maf'
      
      
    # }
    
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
    
    tcga.maf = read.delim(file.path(tcga.maf.dir,tcga.maf.file),stringsAsFactors = F)
    dim(maf) ## 2246  282
    
    print(dim(tcga.maf[tcga.maf$Hugo_Symbol=='TP53',,drop=F]))
    
    table(tcga.maf$Variant_Classification)
    
    tcga.maf.sm = sort(unique(tcga.maf$Tumor_Sample_Barcode))
    
    dim(tcga.maf) ## 2246  282
    tcga.maf = tcga.maf[!is.na(tcga.maf$Variant_Classification) &
                          tcga.maf$Variant_Classification %in% funcs,,drop=F]
    dim(tcga.maf) ## 2238  282
    
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
      load('01_Data/GDAC_data/cnv/TCGA_allCancer.copynumber.all_thresholded_by_gene.txt.RData')
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
  
}

## --------------------------------------------------------

## draw plots 
if(TRUE) {

  ## parse maf
  if(FALSE) {
    
    maf.sampleSummary = getSampleSummary(maf)
    clinical[!clinical$Tumor_Sample_Barcode %in% maf.sampleSummary$Tumor_Sample_Barcode,]
    
    
    plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', 
                   dashboard = TRUE, titvRaw = FALSE)
    
    ## oncoplot for top ten mutated genes.
    oncoplot(maf = maf, top = 20, removeNonMutated = F,
             clinicalFeatures = c('Stage.at.diagnosis'))
    
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
      geom_OneDrive - University of Pittsburghplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
      geom_jitter(shape = 21, size = 5, width = 0.1, height = 0, fill = '#999999') +
      # scale_fill_manual(values = plot.colors) + 
      # stat_compare_means(method = 'wilcox.test') +
      ggtitle('Total TMB (Non-Synonymous Somatic Mutations (SNVs/indels), NSSMs)') +
      theme_classic()
    
    p1
    
    p2 = ggplot(clinical, aes(Location, total)) +
      geom_OneDrive - University of Pittsburghplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
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
      geom_OneDrive - University of Pittsburghplot(width = 0.5, color = '#C0C0C0', size = 3, outlier.shape = NA) +
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
  if(FALSE) {
    
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
      geom_OneDrive - University of Pittsburghplot(width = 0.5, color = '#C0C0C0', size = 2) +
      geom_jitter(width = 0.1, height = 0, size = 2) +
      theme_pubr() +
      stat_compare_means(method = 'wilcox.test')
    
    p1  
    
    
    p1 = ggplot(data.plot, aes(group, t_alt_count)) +
      geom_OneDrive - University of Pittsburghplot(width = 0.5, color = '#C0C0C0', size = 2) +
      geom_jitter(width = 0.1, height = 0, size = 2) +
      theme_pubr() +
      stat_compare_means(method = 'wilcox.test')
    
    p1  
    
    
  }

  ## --------------------------------------------------------
  
  ## DEG analysis or gene set comp: ARID1A mut vs nonmut
  ## DEG analysis: TP53 mut vs nonmut on TCGA
  if(TRUE) {
    
    my.gene = 'TP53' ## ARID1A TP53
    
    if(my.gene == 'ARID1A') {
      clinical$Group.ARID1A = 'nonmut'
      clinical$Group.ARID1A[which(clinical$Tumor_Sample_Barcode %in%
                                    maf$Tumor_Sample_Barcode[
                                      maf$Hugo_Symbol ==my.gene &
                                        maf$Variant_Classification %in% funcs
                                      ])] = 'mut'
      table(clinical$Group.ARID1A)
      # ARID1A.mut ARID1A.nonmut 
      # 4            21
      
      write.csv(clinical,
                file = paste0(clinical.file,'.ARID1A.csv'))
    } else if (my.gene == 'TP53') {
      
      rnaseq.dir = 'results/rnaseq/gene_expr'
      
      atm_expr_norm = read.delim(file.path(rnaseq.dir,'ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.TMM.logCPM.txt'), stringsAsFactors = F, row.names = 1)
      atm_expr_norm_batch = read.delim(file.path(rnaseq.dir,'ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.limmavoomweighted.batchcorr.txt'),stringsAsFactors = F, row.names = 1)
      
      tcga.clinical = data.frame(Tumor_Sample_Barcode = 
                                   tcga.maf.sm,
                                 stringsAsFactors = F)
      dim(tcga.clinical) ## 176
      
      tcga.clinical$Group.TP53 = 'nonmut'
      tcga.clinical$Group.TP53[which(tcga.clinical$Tumor_Sample_Barcode %in%
                                       tcga.maf$Tumor_Sample_Barcode[
                                         tcga.maf$Hugo_Symbol ==my.gene &
                                           tcga.maf$Variant_Classification %in% funcs
                                      ])] = 'mut'
      table(tcga.clinical$Group.TP53)
      # mut nonmut 
      # 118     66
      
      write.csv(tcga.clinical,
                file = paste0(rnaseq.dir,'/TCGA_PAAD.Group.TP53.csv'))
      
      
    }
    
    
    ## --------------------------------------------------------
    
    ## run limma 
    if(TRUE) {
      
      sample.ex = rnaseq.sample.ex
      sample.in = c()
      
      ## settings
      if(TRUE) {
        cancer = "TCGA_PAAD"  ## ATM_PAAD TCGA_PAAD
        fcs = c(2.0,1.5)
        fdrs = c(0.05,0.10)
        rawps = c()
        if(cancer == 'ATM_PAAD') { rawps = c(0.005, 0.01)  }
        limmatypes = c("limmavoom", "limmavoomweighted") ## raw counts (not limma itself!)
        use.spline.flag = 0 ## results are wierd if using splines; NOT using it!
        
        ## ------------------------------------------------------
        
        sample.order = c()
        type = 'coding' ## coding; lincrna
        out.dir = 'diff_expr'
        
        ## ------------------------------------------------------
        
        cpm.flt.flag = 1 ## whether to filter low expr genes by CPM
        cpm.flt.type = 'hard' ## hard or soft
        cpm.string = ''
        if(cpm.flt.type == 'hard') { cpm.thres = 3 } ## 1 2 or 3
        
        # if(cpm.flt.flag == 1) { cpm.string = '.cpm1' }
        
        ## ------------------------------------------------------
        
        # show_col(pal_jco()(10))
        jco.colors = c("#0073C2FF","#EFC000FF","#868686FF","#CD534CFF",
                       "#7AA6DCFF","#003C67FF","#8F7700FF","#3B3B3BFF",
                       "#A73030FF","#4A6990FF")
        
        plot.colors = c('hot' = '#CC0000', 'cold' = '#0000CC', 'med' = '#C0C0C0',
                        'Unknown' = '#000000',
                        'PD' = '#0000CC','SD' = '#00CC00','PR' = '#FF00FF','CR' = '#CC0000',
                        'NA' = '#000000', 'R' = '#CC0000', 'NR' = '#0000CC',
                        # 'A091201_metastatic' = '#EFC000FF','TCGA_primary' = '#0073C2FF',
                        # 'A091201_metastatic' = '#000000','TCGA_primary' = '#C0C0C0',
                        'A091201_metastatic' = '#EFC000FF','TCGA_primary' = '#66CCFF',
                        'Alive' = '#FF0000', 'Dead' = '#000000',
                        'Y' = '#CD534CFF', 'N' = '#7AA6DCFF',
                        'Samples.oneyear.Y' = '#CD534CFF', 'Samples.oneyear.N' = '#7AA6DCFF',
                        'OverOneYear' = '#CD534CFF', 'LessOneYear' = '#7AA6DCFF',
                        'high' = '#000000', 'low' = '#C0C0C0',
                        'C1' = '#E0E0E0', 'C2' = '#A0A0A0', 'C3' = '#606060', 'C4' = '#000000',
                        'ABCS' = '#0073C2FF', 'D' = '#EFC000FF',
                        'ABC' = '#0073C2FF','BC' = '#A73030FF',
                        'A' = '#7AA6DCFF', 'B' = '#003C67FF', 'C' = '#8F7700FF', 'S' = '#3B3B3BFF',
                        'mut' = '#000000', 'nonmut' = '#C0C0C0')
      }

      ## ------------------------------------------------------
      
      ## paths
      path = paste0('results/rnaseq/gene_expr')
      setwd(path)
      
      ## ------------------------------------------------------
      
      ## files
      if(TRUE) {
        
        if(cancer == 'ATM_PAAD') {
          group.file = paste0('../',clinical.file,'.ARID1A.csv')
          expr.file =  'ATM_PAAD.rnaseq.kallisto.raw.txi.coding.txt'
          anno.file = 'gencode.v32.primary_assembly.annotation.maskPAR.gtf.geneinfo'
        } else if (cancer =='TCGA_PAAD') {
          group.file = paste0('TCGA_PAAD.Group.TP53.csv')
          expr.file =  'ATM_PAAD.rnaseq.kallisto.raw.txi.coding.txt'
          anno.file = 'gencode.v32.primary_assembly.annotation.maskPAR.gtf.geneinfo'
          expr.file = '01_Data/GDC_data/rnaseq/rnaseq_counts.TCGA-PAAD.csv'
          anno.file = '01_Data/GDC_data/gencode.gene.info.v22.tsv'
        }
        
      }
      
      ## ------------------------------------------------------
      
      ## import and preprocess data
      if(TRUE) {
        
        if(cancer =='ATM_PAAD') {
          expr = read.delim(expr.file, header = T, stringsAsFactors = F, 
                            row.names = 1)
        } else if(cancer =='TCGA_PAAD') {
          expr = read.csv(expr.file, header = T, stringsAsFactors = F, 
                          row.names = 1)
        }
        
        group = read.csv(group.file, header = T, stringsAsFactors = F, 
                         na.strings = c('na','NA'))
        anno = read.delim(anno.file, header = T, stringsAsFactors = F)
        
        # save(expr, file = paste0(expr.file,'.expr.RData'))
        # load(paste0(expr.file,'.expr.RData'))
        
        if(cancer =='TCGA_PAAD') {
          group$SampleID.RNAseq = gsub('-\\w+{3}-\\w+{4}-\\w+{2}$','',
                                       group$Tumor_Sample_Barcode)
          colnames(expr) = gsub('-\\w+{3}-\\w+{4}-\\w+{2}$','',
                                gsub('[.]','-', colnames(expr)))
        }
        
        group$Group = group[,paste0('Group.',my.gene),drop=F]
        group$Sample = group$SampleID.RNAseq
        group = group[!is.na(group$Sample),]
        
        ## ---------------------------------------------
        
        dim(expr) ## 19784    40
        dim(group) ## ## 40 107
        dim(anno) ## 58395     6
        
        table(group$Group)
        
        anno$gene_id = gsub('[.]\\d+$','',anno$gene_id)
        row.names(anno) = anno$gene_id
        
        expr = merge(anno[,c('gene_name','gene_id')], expr, by = 'row.names')
        row.names(expr) = paste0(expr$gene_name,'!',expr$gene_id)
        expr = expr[,-c(1:3)]
        dim(expr)
        
        ## ---------------------------------------------
        
        ## round numbers to integers (limma will not complain but to meet TMM norm assumption!!)
        dim(expr)
        expr = apply(expr, 2, function(x) round(x, digits = 0))
        dim(expr)
        
        ## ---------------------------------------------
        
        ## remove samples 
        if(length(sample.ex) > 0) {
          expr = expr[,! colnames(expr) %in% sample.ex]
        } else if (length(sample.in) > 0) {
          expr = expr[,colnames(expr) %in% sample.in]
        }
        
        dim(expr) ## 19784    40
        
        row.names(group) = group$Sample
        
        if(length(sample.ex) > 0) {
          group = group[! group$Sample %in% sample.ex,]
        } else if (length(sample.in) > 0) {
          group = group[group$Sample %in% sample.in,]
        }
        
        dim(group) ## 40 107
        
        expr = expr[,colnames(expr) %in% group$Sample]
        group = group[group$Sample %in% colnames(expr),]
        
        ## design matrix colname does not recognize dash or starting with number
        # group$Subject = paste0('P',gsub('-','_',group$Subject)) 
        
        # anno$gene_id = gsub('[.]\\d+$','', anno$gene_id)
        row.names(anno) = anno$gene_id
        
        ## ------------------------------------------------------------------------
        ## Split by gene type 
        ## ------------------------------------------------------------------------
        
        ## split into coding and noncoding
        
        head(anno)
        data.frame(table(anno$gene_type))
        
        anno$gene_id = gsub('[.]\\d+$','', anno$gene_id)
        
        gene.coding = anno[anno$gene_type == 'protein_coding',]
        gene.lincrna = anno[anno$gene_type == 'lincRNA',]
        
        dim(gene.coding)
        dim(gene.lincrna)
        
        expr.coding = expr[row.names(expr) %in% paste0(gene.coding$gene_name,
                                                       '!',gene.coding$gene_id),]
        expr.lincrna = expr[row.names(expr) %in% paste0(gene.lincrna$gene_name,
                                                        '!',gene.lincrna$gene_id),]
        
        write.table(expr.coding,
                    file = paste0(gsub('[.]txt$','',expr.file),
                                  '.',ncol(expr),'sm.coding.matrix'),
                    col.names = NA, row.names = T, sep = '\t', quote = F)
        
        write.table(expr.lincrna,
                    file = paste0(gsub('[.]txt$','',expr.file),
                                  '.',ncol(expr),'sm.lincrna.matrix'),
                    col.names = NA, row.names = T, sep = '\t', quote = F)
        
        ## ---------------------------------------------
        
        expr.original = expr
        
        if(type == 'coding') {
          expr = expr.coding
          
        } else if (type == 'lincrna') {
          expr = expr.lincrna
        }
        
        dim(expr) ## 19784    40
        
        ## ------------------------------------------------------------------------
        ## Select cancer
        ## ------------------------------------------------------------------------
        
        ## sort expr matrix same sample order as group
        expr = expr[,group$Sample]
        print(sum(colnames(expr) == group$Sample) == ncol(expr))
        
        data_sub = expr
        group_sub = group
      }
      
      ## ------------------------------------------------------
      
      ## set up design matrix and contrast
      if(TRUE) {
        
        ## sort matrix sample col to be consistent with group list
        data_sub = data_sub[,group_sub$Sample] 
        # data_sub_print = data.frame(Gene = row.names(data_sub), data_sub)
        # save(data_sub, file=paste0(cancer,".",matrixfile,".data.RData"))
        # write.table(data_sub_print, file=paste0(cancer,".",matrixfile,".data.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
        
        ## filter 
        dim(data_sub)
        
        if(cpm.flt.flag == 0) {
          keep = rowSums(data_sub>0) >= 3 ## at least 3 samples have >0 counts
          data_sub = data_sub[keep,] 
          dim(data_sub)
          ## 18381
          
          save(data_sub, file=paste0(expr.file,'.',ncol(data_sub),'sm.',type,".flt.RData"))
          
          write.table(data_sub, 
                      file=paste0(expr.file,'.',ncol(data_sub),'sm.',type,".flt"), 
                      row.names = T, col.names = NA, quote = F, sep = "\t")
        }
        
        ## ---------------------------------------------
        
        group.col = 'Group'
        out.dir = file.path(out.dir,paste0('Group.',my.gene))
        
        ## ---------------------------------------------
        
        ## prepare design matrix and contrasts
        table(group_sub[,group.col])
        # mut nonmut 
        # 6     16
        
        Group = as.factor(as.vector(unlist(group_sub[,group.col])))
        
        min.sm = NA
        min.sm = min(table(Group))
        if(!is.na(min.sm) & min.sm < 6) { min.sm = 6 }
        # if(!is.na(min.sm) & min.sm < 6) { min.sm = 12 }
        # if(!is.na(min.sm) & min.sm < 12) { min.sm = 12 } ## 02/18/2019: to set min.sm consistent through run=1/2/3....
        print(paste0('Minimum sample per group = ', min.sm))
        
        design = model.matrix(~ 0 + Group)
        
        colnames(design) = gsub('Group','',gsub('Donor','',gsub('RunDate','',colnames(design))))
        head(design)
        
        contrast_matrix = makeContrasts(mut  - nonmut,
                                        levels=design)
        
        if(! dir.exists(out.dir)) { dir.create(out.dir, recursive = T) }
      }
      
      ## ---------------------------------------------
      
      for(limmatype in limmatypes) {
        
        print(paste0("limmatype = ", limmatype))
        output = paste0(out.dir, '/', limmatype, '/', cancer,'.',ncol(data_sub),'sm.',type)
        print(output)
        
        if(! dir.exists(paste0(out.dir, '/', limmatype))) { 
          dir.create(paste0(out.dir, '/', limmatype), recursive = T) 
        }
        
        if( limmatype == "limma" ) { 
          ## create an ExpressionSet obj from a matrix (use Biobase)
          eset = new("ExpressionSet", exprs=as.matrix(data_sub))
          fit = lmFit(eset, design) 
          
        } else if ( limmatype == "limmavoom" | limmatype == "limmavoomweighted") {
          
          dge = DGEList(counts=as.matrix(data_sub))
          dim(dge) ## 19818
          
          ## calculate normalized CPM
          if(11==11) {
            
            dge.norm = calcNormFactors(dge, method="TMM")
            
            write.table(cpm(dge.norm, log=F, prior.count=0.5),
                        file = paste0(out.dir, '/', limmatype, '/',
                                      cancer,'.',ncol(data_sub),'sm.',type,'.TMM.CPM.txt'),
                        sep = '\t', col.names = NA, row.names = T, quote = F)
            
            write.table(cpm(dge.norm, log=T, prior.count=0.5),
                        file = paste0(out.dir, '/', limmatype, '/',
                                      cancer,'.',ncol(data_sub),'sm.',type,'.TMM.logCPM.txt'),
                        sep = '\t', col.names = NA, row.names = T, quote = F)
            
          }
          
          ## do not norm before filtering low expr genes; norm AFTER!
          # dge = calcNormFactors(dge, method="TMM") 
          if(cpm.flt.flag == 1) {
            
            print(paste0('CPM flt type is ',cpm.flt.type, '!'))
            
            ## use cpm > 1 (hard cutoff!)
            if (cpm.flt.type == 'hard') {
              
              dge.cpm = cpm(dge, log=F, prior.count=0.5)
              dge.logcpm = cpm(dge, log=T, prior.count=0.5)
              
              dge.cpm.stats = data.frame(total = 
                                           apply(dge.cpm, 1, 
                                                 function(x) sum(x>cpm.thres, na.rm = T)))
              dge.logcpm.stats = data.frame(total = 
                                              apply(dge.logcpm, 1, 
                                                    function(x) sum(x>-5, na.rm = T)))
              
              ## very stringent low expr filter ... 
              ## requiring all samples having CPM >= 1 - CD8A is lost?!
              # keep = which(dge.cpm.stats$total >= ncol(dge))
              keep = which(dge.cpm.stats$total >= as.integer(min.sm / 2)) 
              
              print(length(keep))
              ## os by 1 year: min.sm = 12; >= 6 samples
              ## os by median: min.sm = 9; >= 4 samples
              ## 02/18/2019: should set all min.sm to 12!!!!!!!!!!!!!!
              
            }
            
            if(cpm.flt.type == 'soft') {
              
              ## https://f1000research.com/articles/5-1408/v1
              # The following filtering rule attempts to keep the maximum number of 
              # interesting genes in the analysis, 
              # but other sensible filtering criteria are also possible. 
              # For example keep <-rowSums(y$counts) > 50 is 
              # a very simple criterion that would keep genes with a total read count of 
              # more than 50. 
              # This would give similar downstream results for this dataset to 
              # the filtering actually used. 
              # Whatever the filtering rule, it should be independent of 
              # the information in the targets file. 
              # It should not make any reference to which RNA libraries belong to 
              # which group, because doing so would bias the subsequent differential expression
              # analysis.
              
              y = dge
              table(rowSums(y$counts == 0) == ncol(y))
              y = y[!rowSums(y$counts == 0) == ncol(y),]
              # at least 2 samples(2 replicates) with a CPM of 10/L, 
              # where L is the # of millions of counts in the smallest library 
              cpm.thres = 10 / min(colSums(y$counts)) * 10^6
              keep = rowSums(cpm(y) > cpm.thres) >= as.integer(min.sm / 2)
              
              dge = y
            }
            
            print(cpm.thres)
            cpm.string = paste0('.cpm', round(cpm.thres, 1))
            output = paste0(output, cpm.string)
            print(output)
            
            ## plot CPM density with threshold
            if(12==12) {
              data.plot = cpm(dge, log=T, prior.count=0.5)
              
              ## plot density, label cutoff line
              data.plot = melt(data.plot)
              colnames(data.plot) = c('Gene', 'Sample', 'logCPM')
              
              plot.title = paste0(expr.file, ' logCPM distribution')
              plot.sub.title = paste0(length(unique(data.plot$Gene)), 
                                      ' genes, ', length(unique(data.plot$Sample)), ' samples')
              p1 = ggplot(data.plot, aes(logCPM)) +
                geom_density(aes(group = Sample), linetype = 'dotted') +
                ggtitle(label = plot.title, subtitle = plot.sub.title) +
                theme_bw() + 
                theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
                theme(axis.text=element_text(size=14), axis.title=element_text(size=12),
                      axis.text.x = element_text(angle = 0),
                      plot.title = element_text(size=12,face="bold", hjust = 0.5)) +  
                theme(strip.text.x = element_text(size =14, face="bold")) +
                geom_vline(xintercept = log(cpm.thres), color = '#CC0000', size = 1)
              
              pdf(file = paste0(output, '.logCPM.density.pdf'), height = 5, width = 5)
              print(p1)
              dev.off()
              
            }
            
            ## recalculate lib size after filtering out low expr genes (Mathew Stephens)
            dim(dge) ## 19818
            dge = dge[keep,,keep.lib.sizes=FALSE]
            dim(dge) ## 16538
            
          } else {
            
            output = paste0(output, cpm.string)
            print(output)
          }
          
          
          ## norm AFTER filtering genes and before voom
          dge = calcNormFactors(dge, method="TMM") 
          
          write.table(cpm(dge, log=F, prior.count=0.5),
                      file = paste0(output,'.TMM.CPM.txt'),
                      sep = '\t', col.names = NA, row.names = T, quote = F)
          
          write.table(cpm(dge, log=T, prior.count=3),
                      file = paste0(output,'.TMM.logCPM.txt'),
                      sep = '\t', col.names = NA, row.names = T, quote = F)
          
          output = paste0(output, ".",limmatype)
          
          
          pdf(file=paste0(output, ".meanvar.pdf"), width=8, height=6)
          
          if( limmatype == "limmavoom") {
            
            # data_voom = voom(dge, normalize.method="none", plot=TRUE) 
            data_voom = voom(dge, design, normalize.method="none", plot=TRUE)
          } else if( limmatype == "limmavoomweighted") {
            data_voom = voomWithQualityWeights(dge, design, normalize.method="none", plot=TRUE)
          }
          dev.off()
          
          save(data_voom, file=paste0(output,".RData"))
          
          ## voom norm matrix
          data_voom_matrix = data.frame(Gene = row.names(data_voom$E), 
                                        round(data_voom$E,2)) 
          
          write.table(data_voom_matrix, file=paste0(output,".txt"), 
                      row.names = F, col.names = T, quote = F, sep = "\t")
          
          
          fit = lmFit(data_voom, design)
          
          
        }
        
        fit2 = contrasts.fit(fit, contrast_matrix)
        fit2 = eBayes(fit2)
        
        coef_colnames = gsub(" - ", "vs", colnames(fit2$cov.coefficients))
        deg_lists = list()
        
        output.0 = output
        
        for (i in 1:length(coef_colnames)) {
          print(paste0("cov.coefficients col number = ", i))
          comp = coef_colnames[i]
          output = paste0(output.0,'.',comp, ".txt")
          table = topTable(fit2, coef=i, adjust="fdr", number=100000, 
                           sort.by="p", resort.by="p", p.value=1, lfc=log2(1.0))
          
          ## add gene name    
          table = data.frame(Gene = row.names(table), table)
          
          ## anti-log FC
          table$FC = 0 
          row_pos = which(table$logFC >= 0)
          row_neg = which(table$logFC < 0)
          table$FC[row_pos] = 2^table$logFC[row_pos]
          table$FC[row_neg] = -2^((-1) * table$logFC[row_neg])
          print(sum(table$FC == 0) == 0)
          head(table)
          write.table(table, output, col.names=T, row.names=F, sep="\t", quote=F)
          
          ## flt DEGs  : FDR
          for(fdr in fdrs) {
            for(fc in fcs) {
              print(paste0(fc, ",", fdr))
              flt.string = paste0(".flt.fdr", fdr, "_fc", fc)
              
              table_flt = table[table$adj.P.Val < fdr & (table$logFC >= log2(fc) | 
                                                           table$logFC <= -log2(fc)),]
              
              if(nrow(table_flt) > 0) { 
                write.table(table_flt, paste0(output,flt.string), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                
                ## print gene symbol only as input for IPA 
                genes = data.frame(Gene = gsub("!\\w+", "", table_flt$Gene), 
                                   table_flt[,colnames(table_flt) %in% 
                                               c("logFC", "P.Value", "adj.P.Val", "FC")])
                print(sum(genes$FC == 0) == 0)
                head(genes)
                genes = genes[genes$Gene != "---",]
                
                write.table(genes, paste0(output, flt.string, ".input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                write.table(genes[genes$logFC > 0, ], paste0(output, flt.string,".up.input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                write.table(genes[genes$logFC < 0, ], paste0(output, flt.string,".down.input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                
                ## --------------------------------------------------
                
                ## heatmap
                genes.sig = table_flt$Gene
                data.plot = data_voom_matrix[data_voom_matrix$Gene%in%genes.sig,]
                row.names(data.plot) = data.plot$Gene
                dim(data.plot)
                
                if(2==2) {
                  
                  ## add annotation 
                  
                  # sample.anno = group_sub[,c('Sample',group.col)]
                  sample.anno = data.frame(Sample = group_sub$Sample,
                                           Group = as.vector(unlist(group_sub[,group.col])))
                    
                  sample.anno[is.na(sample.anno)] = 'NA'
                  
                  ## sort sample anno same as expression matrix
                  sample.anno = sample.anno[order(match(sample.anno$Sample,
                                                        colnames(data.plot))),]
                  row.names(sample.anno) = sample.anno$Sample
                  # sample.anno = sample.anno[,-1]
                  
                  sample.anno.colors = list(
                    Group = plot.colors
                  )
                  
                  plot.anno = HeatmapAnnotation(df = sample.anno[,! colnames(sample.anno)%in%
                                                                   c('Sample'),drop=F],
                                                col = sample.anno.colors)
                  
                }
                
                if(3==3) {
                  
                  data.plot = data.plot[,-1]
                  col.before = colnames(data.plot)
                  col.after = sample.anno[order(match(sample.anno$Sample,
                                                      colnames(data.plot))),]
                  colnames(data.plot) = paste0(col.after$Subject,'|',col.after$Sample)
                  colnames(data.plot) = gsub('_','',gsub('CM_','',colnames(data.plot)))
                  
                  # centered_data = t(scale(t(data.plot), scale=F))
                  centered_data = t(scale(t(data.plot), scale=T))
                  
                  col.title = paste0(cancer,'.', limmatype, '.',comp, 
                                     '.',ncol(data.plot),'sm')
                  row.title = paste0(nrow(data.plot),' Genes (FDR<',fdr, ', FC',fc,')')
                  
                  p2 =Heatmap(centered_data,
                              name = 'log2Expression',
                              column_title = col.title, 
                              row_title = row.title,
                              # column_title_side = 'bottom',
                              column_dend_height = unit(2, "cm"), 
                              row_dend_width = unit(4, "cm"), 
                              # km = 2,
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "ward.D2",
                              clustering_distance_columns = "euclidean",
                              clustering_method_columns = "ward.D2",
                              show_row_names = T,
                              column_title_gp =  gpar(fontsize = 12),
                              bottom_annotation = plot.anno)
                  
                  pdf(file = paste0(output, flt.string,'.heatmap.pdf'), 
                      width = 12, height = 8)
                  print(p2)
                  dev.off()
                  
                }
                
                ## --------------------------------------------------
                
                ## volcano plot 
                data.plot = table[,c('Gene', 'logFC', 'P.Value')]
                data.plot$pass = 'N'
                data.plot$pass[data.plot$P.Value<0.001 & 
                                 (data.plot$logFC>=log2(fc) | data.plot$logFC<=-log2(fc))]='M'
                data.plot$pass[data.plot$Gene %in% 
                                 table_flt$Gene]='Y'
                row.names(data.plot) = data.plot$Gene
                data.plot = data.plot[,-1]
                
                data.plot$Gene = gsub('[!]\\S+$', '', row.names(data.plot))
                
                plot.x.max = round(max(abs(data.plot$logFC)))
                
                if(4==4) {
                  
                  scatter.colors = c('N' = '#C0C0C0', 'Y' = '#CC0000', 'M' = '#CCCC00')
                  scatter.sizes = c('N' = 0.5, 'Y' = 1.5, 'M' = 1.0)
                  
                  plot.title = paste0(cancer,'.', limmatype, '.',comp, 
                                      '.',ncol(data.plot),'sm')
                  plot.sub.title = paste0('Colored in red: ',nrow(table_flt),
                                          ' Genes (FDR<',fdr, ', FC',fc,')')
                  
                  
                  p3 = ggplot(data.plot, aes(logFC, -log10(P.Value))) +
                    geom_point(aes(color = pass, size = pass)) +
                    scale_color_manual(values = scatter.colors) +
                    scale_size_manual(values = scatter.sizes) +
                    # geom_smooth(method = 'lm', se = FALSE, color = 'black', size=1) +
                    xlab('Fold Change (log2)') +
                    ylab('Raw Pvalue (-log10)') +
                    ggtitle(label = plot.title, subtitle = plot.sub.title) +
                    theme_bw() + 
                    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
                    theme(axis.text=element_text(size=14), axis.title=element_text(size=12),
                          plot.title = element_text(size=12,face="bold", hjust = 0.5)) +  
                    theme(strip.text.x = element_text(size =14, face="bold")) +
                    geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = 'dashed', 
                               size = 0.5)+
                    geom_hline(yintercept = -log10(0.001), linetype = 'dashed',
                               size = 0.5) +
                    scale_x_continuous(breaks = seq(-plot.x.max,plot.x.max, by =1),
                                       lim = c(-plot.x.max,plot.x.max))
                  
                  p4 = p3 + 
                    geom_text_repel(
                      data = subset(data.plot, pass == 'Y'),
                      aes(label = Gene),
                      size = 2,
                      box.padding = unit(0.1, "lines"),
                      point.padding = unit(0.1, "lines")
                    )
                  
                  
                  pdf(file = paste0(output, flt.string,'.volcano.pdf'), 
                      width = 6.5, height = 6)
                  print(p3)
                  dev.off()
                  
                  # pdf(file = paste0(output, flt.string,'.volcano.wGeneLabel.pdf'), 
                  #     width = 6.5, height = 6)
                  # print(p4)
                  # dev.off()
                  
                }
                
              }
            }
          }
          
          ## flt DEGs : raw pvalue
          for(rawp in rawps) {
            for(fc in fcs) {
              print(paste0(fc, ",", rawp))
              flt.string = paste0(".flt.rawp", rawp, "_fc", fc)
              
              table_flt = table[table$P.Value < rawp & (table$logFC >= log2(fc) | 
                                                          table$logFC <= -log2(fc)),]
              
              if(nrow(table_flt) > 0) { 
                write.table(table_flt, paste0(output, flt.string), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                
                ## print gene symbol only as input for IPA 
                genes = data.frame(Gene = gsub("!\\w+", "", table_flt$Gene), 
                                   table_flt[,colnames(table_flt) %in% 
                                               c("logFC", "P.Value", "adj.P.Val", "FC")])
                print(sum(genes$FC == 0) == 0)
                head(genes)
                genes = genes[genes$Gene != "---",]
                
                write.table(genes, paste0(output, flt.string, ".input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                write.table(genes[genes$logFC > 0, ], paste0(output, flt.string, ".up.input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                write.table(genes[genes$logFC < 0, ], paste0(output, flt.string,".down.input"), 
                            col.names=T, row.names=F, sep="\t", quote=F)
                
                ## --------------------------------------------------
                
                ## heatmap
                genes.sig = table_flt$Gene
                data.plot = data_voom_matrix[data_voom_matrix$Gene%in%genes.sig,]
                row.names(data.plot) = data.plot$Gene
                dim(data.plot)
                
                
                if(2==2) {
                  
                  ## add annotation 
                  
                  sample.anno = group_sub[,c('Sample',group.col)]
                  sample.anno[is.na(sample.anno)] = 'NA'
                  
                  ## sort sample anno same as expression matrix
                  sample.anno = sample.anno[order(match(sample.anno$Sample,
                                                        colnames(data.plot))),]
                  row.names(sample.anno) = sample.anno$Sample
                  # sample.anno = sample.anno[,-1]
                  
                    sample.anno.colors = list(
                      Group = plot.colors
                    )

                  
                  plot.anno = HeatmapAnnotation(df = sample.anno[,! colnames(sample.anno)%in%
                                                                   c('Sample'),drop=F],
                                                col = sample.anno.colors)
                  
                }
                
                if(3==3) {
                  
                  data.plot = data.plot[,-1]
                  
                  # centered_data = t(scale(t(data.plot), scale=F))
                  centered_data = t(scale(t(data.plot), scale=T))
                  
                  col.title = paste0(cancer,'.', limmatype, '.',comp, 
                                     '.',ncol(data.plot),'sm')
                  row.title = paste0(nrow(data.plot),' Genes (rawP<',rawp, ', FC',fc,')')
                  
                  p2 =Heatmap(centered_data,
                              name = 'log2Expression',
                              column_title = col.title, 
                              row_title = row.title,
                              # column_title_side = 'bottom',
                              column_dend_height = unit(2, "cm"), 
                              row_dend_width = unit(4, "cm"), 
                              # km = 2,
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "ward.D2",
                              clustering_distance_columns = "euclidean",
                              clustering_method_columns = "ward.D2",
                              show_row_names = T,
                              column_title_gp =  gpar(fontsize = 12),
                              bottom_annotation = plot.anno)
                  
                  pdf(file = paste0(output, flt.string,'.heatmap.pdf'), 
                      width = 12, height = 8)
                  print(p2)
                  dev.off()
                  
                }
                
                ## --------------------------------------------------
                
                ## volcano plot 
                data.plot = table[,c('Gene', 'logFC', 'P.Value')]
                data.plot$pass = 'N'
                data.plot$pass[data.plot$Gene %in% 
                                 table_flt$Gene]='Y'
                row.names(data.plot) = data.plot$Gene
                data.plot = data.plot[,-1]
                
                data.plot$Gene = gsub('[!]\\S+$', '', row.names(data.plot))
                
                plot.x.max = round(max(abs(data.plot$logFC)))
                
                if(4==4) {
                  
                  scatter.colors = c('N' = '#C0C0C0', 'Y' = '#CC0000')
                  scatter.sizes = c('N' = 0.5, 'Y' = 1.5)
                  
                  plot.title = paste0(cancer,'.', limmatype, '.',comp, 
                                      '.',ncol(data.plot),'sm')
                  plot.sub.title = paste0('Colored in red: ',nrow(table_flt),
                                          ' Genes (rawP<',rawp, ', FC',fc,')')
                  
                  
                  p3 = ggplot(data.plot, aes(logFC, -log10(P.Value))) +
                    geom_point(aes(color = pass, size = pass)) +
                    scale_color_manual(values = scatter.colors) +
                    scale_size_manual(values = scatter.sizes) +
                    # geom_smooth(method = 'lm', se = FALSE, color = 'black', size=1) +
                    xlab('Fold Change (log2)') +
                    ylab('Raw Pvalue (-log10)') +
                    ggtitle(label = plot.title, subtitle = plot.sub.title) +
                    theme_bw() + 
                    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
                    theme(axis.text=element_text(size=14), axis.title=element_text(size=12),
                          plot.title = element_text(size=12,face="bold", hjust = 0.5)) +  
                    theme(strip.text.x = element_text(size =14, face="bold")) +
                    geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = 'dashed', 
                               size = 0.5)+
                    geom_hline(yintercept = -log10(rawp), linetype = 'dashed', 
                               size = 0.5) +
                    scale_x_continuous(breaks = seq(-plot.x.max,plot.x.max, by =1),
                                       lim = c(-plot.x.max,plot.x.max))
                  
                  p4 = p3 + 
                    geom_text_repel(
                      data = subset(data.plot, pass == 'Y'),
                      aes(label = Gene),
                      size = 2,
                      box.padding = unit(0.1, "lines"),
                      point.padding = unit(0.1, "lines")
                    )
                  
                  
                  pdf(file = paste0(output, flt.string,'.volcano.pdf'), 
                      width = 6.5, height = 6)
                  print(p3)
                  dev.off()
                  
                  # pdf(file = paste0(output, flt.string,'.volcano.wGeneLabel.pdf'), 
                  #     width = 6.5, height = 6)
                  # print(p4)
                  # dev.off()
                  
                }
              }
            }
          }
        }
        
        ## pca plot 
        if(1==1) {
          
          data.plot = t(data_voom_matrix[,-1])
          data.plot = merge(data.plot,group_sub[,c(group.col),drop=F], 
                            by = 'row.names')
          row.names(data.plot) = data.plot[,1]
          data.plot = data.plot[,-1]
          
          data.pca = prcomp(data.plot[,1:nrow(data_voom_matrix)])
          
          
            p5 = autoplot(data.pca, data = data.plot, colour = 'Group',
                          label = TRUE, label.size = 3)
          
          
          pdf(file = paste0(output, '.',ncol(data.plot),'genes_',nrow(data.plot),'sm',
                            '.pca.pdf'), 
              width = 6.5, height = 5)
          print(p5)
          dev.off()
          
        }
      }
      
    }
    
    ## --------------------------------------------------------
    
    ## compare hallmark gene set scores 
    if(TRUE) {
      
      ## import hallmark genesets
      if(TRUE) {
        hallmark.file = '../../../h.all.v7.3.symbols.gmt'
        
        hallmark = read.delim(hallmark.file,
                              stringsAsFactors = F, header = F)
        hallmark = data.frame(t(hallmark),stringsAsFactors = F)
        colnames(hallmark) = hallmark[1,]
        
        hallmark = hallmark[-c(1:2),]
        
        rm(x)
        x =list()
        for(i in 1:ncol(hallmark)) {
          
          x[[i]] = as.vector(hallmark[,i])
          x[[i]] = x[[i]][x[[i]] !='']
          
        }
        names(x) = colnames(hallmark)
        hallmark = x
      }
      
      ## --------------------------------------------------------
      
      ## gsva 
      if(TRUE) {
        
        library(GSVA)
        expr.file = 'diff_expr/Group.ARID1A/limmavoomweighted/ATM_PAAD.23sm.coding.cpm3.TMM.logCPM.txt'
        expr = read.delim(expr.file, stringsAsFactors = F)
        expr$X = gsub('[!]\\S+$','',expr$X)
        rm(x)
        x = as.character(expr$X)
        expr = as.matrix(expr[,-1])
        rownames(expr) = x
        
        gsva_es = gsva(expr,
                       hallmark,
                       min.sz = 15,
                       max.sz = 500,
                       mx.diff=T,
                       method = 'gsva'
        )
        
        write.csv(gsva_es,
                  file = paste0(expr.file,'.hallmark.gsva.csv'))
        
        rm(x)
        x = reshape2::melt(gsva_es)
        p1 = ggplot(x, aes(value)) + geom_density(aes(color = Var2))
        p1
        
      }
      
      ## --------------------------------------------------------

      ## compare hallmark gene sets between mut vs nonmut
      if(TRUE) {
        
        rm(x)
        x = clinical
        x$Sample = x$SampleID.RNAseq
        x = x[!is.na(x$Sample),]
        
        data.stats = NULL 
        
        for(i in 1:nrow(gsva_es)) {
          print(my.set)
          my.set = row.names(gsva_es)[i]
          my.df = NULL 
          my.df = t(gsva_es[i,,drop=F])
          
          my.df = data.frame(Sample = row.names(my.df),
                             my.df)
          
          my.df = merge(my.df, x[,c('Sample','Group.ARID1A')], by = 'Sample')
          
          data.stats = rbind(data.stats,
                             tidy(t.test(
                               my.df[my.df$Group.ARID1A=='mut',2],
                               my.df[my.df$Group.ARID1A=='nonmut',2]
                             )))
          
          
        }
        
        data.stats = data.frame(data.stats, stringsAsFactors = F)
        
        data.stats$p.adj = p.adjust(data.stats$p.value,method = 'fdr')
        
        row.names(data.stats) = row.names(gsva_es)
        
        write.csv(data.stats,
                  file = paste0(expr.file,'.hallmark.gsva.ttest.csv'))
        
      }
      
    }
  }
  
  ## --------------------------------------------------------
  
  ## CDH1: spread of sample-wise gene expression value in ATM vs TCGA (set as same unit) 
  if(TRUE) {
    
    ## paths
    path = paste0('results/rnaseq/gene_expr')
    setwd(path)
    
    ## ------------------------------------------------------
    
    ## norm and batch correction between TCGA and ATM cohorts 
    if(TRUE) {
      
      ## files
      if(TRUE) {
        atm.expr.file =  'ATM_PAAD.rnaseq.kallisto.raw.txi.coding.txt'
        atm.anno.file = 'gencode.v32.primary_assembly.annotation.maskPAR.gtf.geneinfo'
        
        tcga.expr.file = '01_Data/GDC_data/rnaseq/rnaseq_counts.TCGA-PAAD.csv'
        tcga.anno.file = '01_Data/GDC_data/gencode.gene.info.v22.tsv'
      }
      
      ## ------------------------------------------------------
      
      ## import data 
      if(TRUE) {
        
        atm.expr = read.delim(atm.expr.file, stringsAsFactors = F,row.names = 1)
        atm.anno = read.delim(atm.anno.file, stringsAsFactors = F)
        tcga.expr = read.csv(tcga.expr.file, stringsAsFactors = F,row.names = 1)
        tcga.anno = read.delim(tcga.anno.file, stringsAsFactors = F)
        
      }
      
      ## ------------------------------------------------------
      
      ## preprocess
      if(TRUE) {
        
        dim(atm.expr)
        row.names(atm.expr) = gsub('^\\S+[!]','',row.names(atm.expr))
        
        dim(tcga.expr)
        
        atm.expr = merge(atm.expr,
                         tcga.expr, by = 'row.names')
        dim(atm.expr) ## 19600   207
        
        colnames(atm.expr)[1] = 'gene_id'
        atm.anno$gene_id = gsub('[.]\\d+$','',atm.anno$gene_id)
        
        atm.expr = merge(atm.anno[,c('gene_id','gene_name')], 
                         atm.expr, by = 'gene_id')
        
        row.names(atm.expr) = paste0(atm.expr$gene_name,'!',atm.expr$gene_id)
        atm.expr = atm.expr[,-c(1:2)]
        
        dim(atm.expr) ## 19600   206
        
        atm.expr = atm.expr[,!colnames(atm.expr) %in% rnaseq.sample.ex]
        dim(atm.expr) ## 19600   205
        
        write.csv(atm.expr,
                  paste0(atm.expr.file,'.wTCGA.csv'))
        
      }
      
      ## ---------------------------------------------------
      
      ## norm and batch correct 
      if(TRUE) {
        
        data_sub = NULL 
        group_sub = NULL 
        
        data_sub = atm.expr 
        rm(x,y)
        x = grep('^TCGA',colnames(atm.expr))
        y  = which(!colnames(atm.expr) %in% colnames(atm.expr)[x])
        group_sub = data.frame(Sample = colnames(atm.expr),
                               Cohort =NA,
                               stringsAsFactors = F) 
        group_sub$Cohort[x] = 'TCGA_PAAD'
        group_sub$Cohort[y]= 'ATM_PAAD'
        group_sub_batch = factor(group_sub$Cohort)
        min.sm = min(table(group_sub$Cohort))
        print(min.sm)
        
        all.equal(group_sub$Sample, colnames(data_sub))
        
        ## ---------------------------------------------------
        
        ## run limma 
        if(TRUE) {
          
          cpm.flt.flag = 1
          cpm.flt.type = 'hard'
          cpm.thres = 3
          
          cancer = 'ATM_PAAD_TCGA_PAAD'
          limmatype = 'limmavoomweighted'
          type = 'coding'
          
          ## ---------------------------------------------------
          
          
          print(paste0("limmatype = ", limmatype))
          output = paste0(cancer,'.',ncol(data_sub),'sm.',type)
          print(output)
          
          ## ---------------------------------------------------
          
          dge = DGEList(counts=as.matrix(data_sub))
          dim(dge) ##  19883   228
          
          ## calculate normalized CPM
          if(TRUE) {
            
            dge.norm = calcNormFactors(dge, method="TMM")
            
            write.table(cpm(dge.norm, log=F, prior.count=3),
                        file = paste0(output,'.TMM.CPM.txt'),
                        sep = '\t', col.names = NA, row.names = T, quote = F)
            
            write.table(cpm(dge.norm, log=T, prior.count=3),
                        file = paste0(output,'.TMM.logCPM.txt'),
                        sep = '\t', col.names = NA, row.names = T, quote = F)
            
          }
          
          ## ---------------------------------------------------
          
          ## do not norm before filtering low expr genes; norm AFTER!
          # dge = calcNormFactors(dge, method="TMM") 
          if(cpm.flt.flag == 1) {
            
            print(paste0('CPM flt type is ',cpm.flt.type, '!'))
            
            ## use cpm > 1 (hard cutoff!)
            if (cpm.flt.type == 'hard') {
              
              dge.cpm = cpm(dge, log=F, prior.count=3)
              dge.logcpm = cpm(dge, log=T, prior.count=3)
              
              dge.cpm.stats = data.frame(total = 
                                           apply(dge.cpm, 1, 
                                                 function(x) sum(x>cpm.thres, na.rm = T)))
              dge.logcpm.stats = data.frame(total = 
                                              apply(dge.logcpm, 1, 
                                                    function(x) sum(x>-5, na.rm = T)))
              
              ## very stringent low expr filter ... 
              ## requiring all samples having CPM >= 1 - CD8A is lost?!
              # keep = which(dge.cpm.stats$total >= ncol(dge))
              keep = which(dge.cpm.stats$total >= as.integer(min.sm / 2)) 
              
              print(length(keep)) ## 15301
              
            }
            

            if(cpm.flt.type == 'soft') {
              
              ## https://f1000research.com/articles/5-1408/v1
              # The following filtering rule attempts to keep the maximum number of 
              # interesting genes in the analysis, 
              # but other sensible filtering criteria are also possible. 
              # For example keep <-rowSums(y$counts) > 50 is 
              # a very simple criterion that would keep genes with a total read count of 
              # more than 50. 
              # This would give similar downstream results for this dataset to 
              # the filtering actually used. 
              # Whatever the filtering rule, it should be independent of 
              # the information in the targets file. 
              # It should not make any reference to which RNA libraries belong to 
              # which group, because doing so would bias the subsequent differential expression
              # analysis.
              
              y = dge
              table(rowSums(y$counts == 0) == ncol(y))
              y = y[!rowSums(y$counts == 0) == ncol(y),]
              # at least 2 samples(2 replicates) with a CPM of 10/L, 
              # where L is the # of millions of counts in the smallest library 
              cpm.thres = 10 / min(colSums(y$counts)) * 10^6
              keep = rowSums(cpm(y) > cpm.thres) >= as.integer(min.sm / 2)
              
              dge = y
            }
            
            ## plot CPM density with threshold
            if(TRUE) {
              data.plot = cpm(dge, log=T, prior.count=3)
              
              ## plot density, label cutoff line
              data.plot = melt(data.plot)
              colnames(data.plot) = c('Gene', 'Sample', 'logCPM')
              
              plot.title = paste0(output, ' logCPM distribution')
              plot.sub.title = paste0(length(unique(data.plot$Gene)), 
                                      ' genes, ', length(unique(data.plot$Sample)), ' samples')
              p1 = ggplot(data.plot, aes(logCPM)) +
                geom_density(aes(group = Sample), linetype = 'dotted') +
                ggtitle(label = plot.title, subtitle = plot.sub.title) +
                theme_bw() + 
                theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
                theme(axis.text=element_text(size=14), axis.title=element_text(size=12),
                      axis.text.x = element_text(angle = 0),
                      plot.title = element_text(size=12,face="bold", hjust = 0.5)) +  
                theme(strip.text.x = element_text(size =14, face="bold")) +
                geom_vline(xintercept = log(cpm.thres), color = '#CC0000', size = 1)
              
              pdf(file = paste0(output, '.logCPM.density.pdf'), height = 5, width = 5)
              print(p1)
              dev.off()
              
            }
            
            ## recalculate lib size after filtering out low expr genes (Mathew Stephens)
            dim(dge) ## 19883
            dge = dge[keep,,keep.lib.sizes=FALSE]
            dim(dge) ## 15301
            
          } else {
            
            # output = paste0(output, cpm.string)
            # print(output)
          }
          
          ## ---------------------------------------------------
          
          print(cpm.thres)
          cpm.string = paste0('.cpm', round(cpm.thres, 1))
          output = paste0(output, cpm.string)
          print(output)
          
          ## ---------------------------------------------------
          
          ## norm AFTER filtering genes and before voom
          dge = calcNormFactors(dge, method="TMM") 
          
          ## produce new CPM and logCPM matrix after filtering
          write.table(cpm(dge, log=F, prior.count=3),
                      file = paste0(output,'.TMM.CPM.txt'),
                      sep = '\t', col.names = NA, row.names = T, quote = F)
          write.table(cpm(dge, log=T, prior.count=3),
                      file = paste0(output,'.TMM.logCPM.txt'),
                      sep = '\t', col.names = NA, row.names = T, quote = F)
          
          output = paste0(output, '.',limmatype)
          print(output)
          
          ## ---------------------------------------------------
          
          pdf(file=paste0(output, ".meanvar.pdf"), width=8, height=6)
          
          if( limmatype == "limmavoom") {
            
            # data_voom = voom(dge, normalize.method="none", plot=TRUE) 
            data_voom = voom(dge, normalize.method="none", plot=TRUE)
          } else if( limmatype == "limmavoomweighted") {
            # data_voom = voomWithQualityWeights(dge, design, normalize.method="none", plot=TRUE)
            data_voom = voomWithQualityWeights(dge, normalize.method="none", plot=TRUE)
          }
          dev.off()
          
          save(data_voom, file=paste0(output,".RData"))
          # save(design, contrast_matrix, data_voom, file=paste0(output,".forLei.RData"))
          
          ## ---------------------------------------------------
          
          ## voom norm matrix
          data_voom_matrix = data.frame(Gene = row.names(data_voom$E), 
                                        round(data_voom$E,2)) 
          
          write.table(data_voom_matrix, file=paste0(output,".txt"), 
                      row.names = F, col.names = T, quote = F, sep = "\t")
          
          ## ---------------------------------------------------
          
          ## correct for batch for heatmaps (NOT for DEG!!!)
          if(TRUE) {
            
            data_voom_batch = removeBatchEffect(data_voom, batch = group_sub_batch)
            dim(data_voom_batch) ## 15301   228
            
            write.table(data_voom_batch, file=paste0(output,".batchcorr.txt"), 
                        row.names = T, col.names = NA, quote = F, sep = "\t")
          }
          
          ## ---------------------------------------------------
          
          # fit = lmFit(data_voom, design)
          #
          # fit2 = contrasts.fit(fit, contrast_matrix)
          # fit2 = eBayes(fit2)
          # 
          # coef_colnames = gsub(" - ", "vs", colnames(fit2$cov.coefficients))
          
          ## ---------------------------------------------------
          
          ## pca plot, before and after batch correction
          if(TRUE) {
            
            row.names(group_sub) = group_sub$Sample
            
            ## before batch correction 
            data.plot = t(data_voom_matrix[,-1])
            data.plot = merge(data.plot, group_sub[,-1,drop=F], by = 'row.names')
            row.names(data.plot) = data.plot[,1]
            data.plot = data.plot[,-1]
            
            data.pca = prcomp(data.plot[,1:nrow(data_voom_matrix)])
            p5 = autoplot(data.pca, data = data.plot, colour = 'Cohort',
                          label = FALSE, label.size = 3,
                          frame = TRUE, frame.type = 'norm',frame.colour = 'Cohort') +
              theme_minimal() +
              ggtitle(paste0(output), subtitle = paste0('before batch correction'))
            
            pdf(file = paste0(output, '.',ncol(data.plot),'genes_',nrow(data.plot),'sm',
                              '.pca.pdf'), 
                width = 6.5, height = 5)
            print(p5)
            dev.off()
            
            ## ---------------------------------------------------
            
            ## after  batch correction 
            data.plot = t(data_voom_batch)
            data.plot = merge(data.plot, group_sub[,-1,drop=F], by = 'row.names')
            row.names(data.plot) = data.plot[,1]
            data.plot = data.plot[,-1]
            
            data.pca = prcomp(data.plot[,1:nrow(data_voom_batch)])
            p6 = autoplot(data.pca, data = data.plot, colour = 'Cohort',
                          label = FALSE, label.size = 3,
                          frame = TRUE, frame.type = 'norm',frame.colour = 'Cohort') +
              theme_minimal() +
              ggtitle(paste0(output), subtitle = paste0('after batch correction'))
            
            pdf(file = paste0(output, '.',ncol(data.plot),'genes_',nrow(data.plot),'sm',
                              '.batchcorr.pca.pdf'), 
                width = 6.5, height = 5)
            print(p6)
            dev.off()
            
          }
          
          
        }
      
      }
      
      
    }

    ## ---------------------------------------------------
    
    atm_expr_norm = read.delim('ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.TMM.logCPM.txt', stringsAsFactors = F, row.names = 1)
    atm_expr_norm_batch = read.delim('ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.limmavoomweighted.batchcorr.txt',stringsAsFactors = F, row.names = 1)
    
    ## ---------------------------------------------------
    
    dim(atm_expr_norm)
    dim(atm_expr_norm_batch)
    
    ## CDH1: spread of sample-wise gene expression value in ATM vs TCGA (set as same unit) 
    if(TRUE) {
      
      data.plot = NULL 
      data.plot = atm_expr_norm ## atm_expr_norm atm_expr_norm_batch
      
      my.genes = c('TP53')
      
      data.plot = data.frame(t(
        data.plot[gsub('[!]\\S+$','',row.names(data.plot)) %in% my.genes,,
                            drop=F]))
      
      colnames(data.plot) = gsub('[.]\\S+$','',colnames(data.plot))
      
      data.plot = data.frame(Sample = row.names(data.plot),
                             data.plot,
                             stringsAsFactors = F)
      
      data.plot$Group = NA
      data.plot$Group[grep('^TCGA',data.plot$Sample)] = 'TCGA_PAAD'
      data.plot$Group[-grep('^TCGA',data.plot$Sample)] = 'ATM_PAAD'
      table(data.plot$Group)
      
      p1 = ggplot(data.plot, aes(Group, TP53)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, height = 0) +
        ggtitle('before batchcorrection') +
        # ggtitle('after batchcorrection') +
        theme_pubr()
      
      p1
      
      ## --------------------------------------------------------
      
      ## randomly pick 50 genes from the expr matrix .....
      
      data.plot = NULL 
      data.plot = atm_expr_norm ## atm_expr_norm atm_expr_norm_batch
      
      my.genes = c('TP53')
      
      rm(x)
      set.seed(10)
      x = sample(1:nrow(data.plot), 50)
      x
      
      data.plot = data.frame(t(
        data.plot[gsub('[!]\\S+$','',row.names(data.plot)) %in% 
                    unique(c(my.genes, gsub('[!]\\S+$','',row.names(data.plot))[x])),,
                  drop=F]))
      
      colnames(data.plot) = gsub('[.]\\S+$','',colnames(data.plot))
      
      data.plot = data.frame(Sample = row.names(data.plot),
                             data.plot,
                             stringsAsFactors = F)
      
      data.plot$Group = NA
      data.plot$Group[grep('^TCGA',data.plot$Sample)] = 'TCGA_PAAD'
      data.plot$Group[-grep('^TCGA',data.plot$Sample)] = 'ATM_PAAD'
      table(data.plot$Group)
      
      data.plot = reshape2::melt(data.plot)
      colnames(data.plot)[3:4] = c('Gene','Expr')
      
      p2 = ggplot(data.plot, aes(Gene, Expr)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, height = 0) +
        ggtitle('before batchcorrection') +
        # ggtitle('after batchcorrection') +
        theme_pubr(x.text.angle = 90) +
        facet_grid(Group~.)
      
      p2
      
    }
    
  }
  
  ## --------------------------------------------------------
  
  ## DEG analysis: P53 target molecules : compare ATM P53 output (nonmutant tumors) vs TCGA P53 non mutant samples vs TCGA P53 mutant samples 

  ## set path back! 
  path = 'results/rnaseq/gene_expr'
  setwd(path)
  
  my.run = 4 ## 1=200 targets, before 2=200 targets, after 3=92 targets, before 4=92 targets, after 
  my.direction = 'down' ## up or down ; only applicable for run 3 and 4
  
  if(TRUE) {
   
    ## prep input
    if(TRUE) {
      
      tp53.genes = read.delim('HALLMARK_P53_PATHWAY.txt',stringsAsFactors = F,skip = 1)
      
      colnames(tp53.genes)[1] = 'HALLMARK_P53_PATHWAY'
      
      ## --------------------------------------------------------
      
      ## whether to further narrow down on genes....
      if(my.run %in% c(3,4)) {
      
        tp53.degs = read.delim('diff_expr/Group.TP53/limmavoom/TCGA_PAAD.176sm.coding.cpm3.limmavoom.mutvsnonmut.txt',stringsAsFactors = F)
        tp53.degs.up = tp53.degs[tp53.degs$logFC > 0 & tp53.degs$adj.P.Val<0.05,]
        tp53.degs.down = tp53.degs[tp53.degs$logFC < 0 &tp53.degs$adj.P.Val<0.05,]
        
        ## --------------------------------------------------------
        
        tp53.genes[tp53.genes$HALLMARK_P53_PATHWAY %in% 
                     gsub('[!]\\S+$','',tp53.degs.up$Gene),]
        
        tp53.genes[tp53.genes$HALLMARK_P53_PATHWAY %in% 
                     gsub('[!]\\S+$','',tp53.degs.down$Gene),]

        ## --------------------------------------------------------
        
        dim(tp53.genes) ## 200
        # tp53.genes = tp53.genes[tp53.genes$HALLMARK_P53_PATHWAY %in% 
        #                           gsub('[!]\\S+$','',tp53.degs.down$Gene),,drop=F]
        tp53.genes$direction = 'unchanged'
        tp53.genes$direction[tp53.genes$HALLMARK_P53_PATHWAY %in% 
                               gsub('[!]\\S+$','',tp53.degs.up$Gene)] = 'up'
        tp53.genes$direction[tp53.genes$HALLMARK_P53_PATHWAY %in% 
                               gsub('[!]\\S+$','',tp53.degs.down$Gene)] = 'down'
        dim(tp53.genes) ## 200
        table(tp53.genes$direction)
        # down unchanged        up 
        # 26       108        66
        
      }
      
      ## --------------------------------------------------------
      
      tp53.genes[tp53.genes=='TP53'] ## TP53
      
      ## --------------------------------------------------------
      
      rnaseq.dir = 'results/rnaseq/gene_expr'
      
      atm_expr_norm = read.delim(file.path(rnaseq.dir,'ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.TMM.logCPM.txt'), stringsAsFactors = F, row.names = 1)
      atm_expr_norm_batch = read.delim(file.path(rnaseq.dir,'ATM_PAAD_TCGA_PAAD.205sm.coding.cpm3.limmavoomweighted.batchcorr.txt'),stringsAsFactors = F, row.names = 1)
      
      tcga.clinical = read.csv(paste0(rnaseq.dir,'/TCGA_PAAD.Group.TP53.csv'), stringsAsFactors = F)
      
    }
    
    ## --------------------------------------------------------
    
    ## format input 
    if(TRUE) {
      data.plot = NULL 
      if(my.run %in% c(1,3)) {
        data.plot = atm_expr_norm[gsub('[!]\\S+$','',row.names(atm_expr_norm_batch)) %in% tp53.genes$HALLMARK_P53_PATHWAY,,drop=F]
      } else if(my.run %in% c(2,4)) {
        data.plot = atm_expr_norm_batch[gsub('[!]\\S+$','',row.names(atm_expr_norm_batch)) %in% tp53.genes$HALLMARK_P53_PATHWAY,,drop=F]
      }
      
      dim(data.plot)
      
      data.plot.0 = NULL 
      data.plot.0 = data.plot
      colnames(data.plot.0) = gsub('[.]\\w{3}[.]\\w{4}[.]\\w{2}$','',
                                   colnames(data.plot.0))
      
      if (my.run %in% c(3,4)) {
        table(tp53.genes$direction)
        rm(x,y)
        x = data.plot[gsub('[!]\\S+$','',row.names(data.plot)) %in%
                        tp53.genes$HALLMARK_P53_PATHWAY[tp53.genes$direction=='up'],,drop=F]
        dim(x) ## 66
        y = data.plot[gsub('[!]\\S+$','',row.names(data.plot)) %in%
                        tp53.genes$HALLMARK_P53_PATHWAY[tp53.genes$direction=='down'],,drop=F]
        dim(y) ## 26
        
        # data.plot = rbind(x,y)
        if(my.direction == 'up') { data.plot = x}
        if(my.direction == 'down') { data.plot = y}
        
      }
      
      data.plot.0 = data.plot.0[row.names(data.plot.0) %in% 
                                  row.names(data.plot),,drop=F]
    }

    ## --------------------------------------------------------
    
    ## Tp53 signaling score 
    if(TRUE) {
      
      data.plot = data.frame(TP53.score = apply(data.plot,2,mean,na.rm=T))
      
      data.plot = data.frame(Sample = row.names(data.plot),
                             data.plot,
                             stringsAsFactors = F)
      
      data.plot$Cohort = NA
      data.plot$Cohort[grep('^TCGA',data.plot$Sample)] = 'TCGA_PAAD'
      data.plot$Cohort[-grep('^TCGA',data.plot$Sample)] = 'ATM_PAAD'
      table(data.plot$Cohort)
      dim(data.plot) ## 205
      
      data.plot$Sample = gsub('-\\w{3}-\\w{4}-\\w{2}$','',
                              gsub('[.]','-',data.plot$Sample))
      
      data.plot = data.plot[data.plot$Sample %in%
                              c(clinical$SampleID.RNAseq,
                                gsub('-\\w{3}-\\w{4}-\\w{2}$','',
                                     tcga.clinical$Tumor_Sample_Barcode)),]
      
      data.plot$Group.TP53 = 'nonmut'
      rm(x,y)
      x = clinical$SampleID.RNAseq[clinical$SampleID.DNAseq %in%
                                     gsub('-tumor','',
                                          maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=='TP53'])]
      y = gsub('-\\w{3}-\\w{4}-\\w{2}$','',tcga.maf$Tumor_Sample_Barcode[tcga.maf$Hugo_Symbol=='TP53'])
      data.plot$Group.TP53[data.plot$Sample %in% c(x,y)] = 'mut'
      table(data.plot$Group.TP53)
      
      data.plot$Cohort_Group = paste0(data.plot$Cohort,'_',data.plot$Group.TP53)
      data.plot = reshape2::melt(data.plot)
      colnames(data.plot)[(ncol(data.plot)-1):ncol(data.plot)] = c('Gene','Expr')
      
      if(my.run %in% c(1,2)) {
        plot.title = paste0('TP53 signaling: ',nrow(tp53.genes),' genes')
      } else if (my.run %in% c(3,4)) {
        plot.title = paste0('TP53 signaling: ',nrow(tp53.genes[
          tp53.genes$direction==my.direction,,drop=F
          ]),' genes ',my.direction)
      }
      
      if(my.run %in% c(1,3)) {
        plot.subtitle = 'before batchcorrection'
      } else if (my.run %in% c(2,4)) {
        plot.subtitle = 'after batchcorrection'
      }
      
      p2 = NULL 
      p2 = ggplot(data.plot, aes(Cohort_Group, Expr)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, height = 0) +
        ggtitle(plot.title, subtitle = plot.subtitle) +
        theme_pubr(x.text.angle = 90) +
        stat_compare_means(method = 't.test',
                           comparisons = list(c('ATM_PAAD_nonmut','TCGA_PAAD_mut'),
                                              c('ATM_PAAD_nonmut','TCGA_PAAD_nonmut'),
                                              c('TCGA_PAAD_mut','TCGA_PAAD_nonmut')))
      print(p2)
      
    }
    
    ## --------------------------------------------------------
    
    dim(data.plot.0)
    table(data.plot$Cohort_Group)
    
    data.plot$Sample = gsub('[-]','.',data.plot$Sample)
    
    data.plot.0 = data.plot.0[,data.plot$Sample]
    
    ## TP53 target heatmap 
    if(TRUE) {
      
      if(2==2) {
        
        ## add annotation 
        sample.anno = data.plot
        
        sample.anno[is.na(sample.anno)] = 'NA'
        
        ## sort sample anno same as expression matrix
        sample.anno = sample.anno[order(match(sample.anno$Sample,
                                              colnames(data.plot.0))),]
        row.names(sample.anno) = sample.anno$Sample
        # sample.anno = sample.anno[,-1]
        
        sample.anno.colors = list(
          Cohort = plot.colors,
          Group.TP53 = plot.colors
        )
        
        plot.anno = HeatmapAnnotation(df = sample.anno[,c('Cohort','Group.TP53'),drop=F],
                                      col = sample.anno.colors)
        
      }
      
      if(3==3) {

        centered_data = NULL 
        centered_data = t(scale(t(data.plot.0), scale=T))
        row.names(centered_data) = gsub('[!]\\S+$','',row.names(centered_data))
        
        col.title = paste0(plot.title,' ',plot.subtitle,
                           '.',ncol(centered_data),'sm')
        row.title = paste0(nrow(centered_data),' Genes')
        
        p3 = NULL 
        p3 =Heatmap(centered_data,
                    name = 'log2Expression',
                    column_title = col.title, 
                    row_title = row.title,
                    # column_title_side = 'bottom',
                    column_dend_height = unit(2, "cm"), 
                    row_dend_width = unit(4, "cm"), 
                    # km = 2,
                    clustering_distance_rows = "euclidean",
                    clustering_method_rows = "ward.D2",
                    clustering_distance_columns = "euclidean",
                    clustering_method_columns = "ward.D2",
                    show_row_names = F,
                    show_column_names = F,
                    column_title_gp =  gpar(fontsize = 12),
                    top_annotation = plot.anno)

        print(p3)
        
      }
    }
    
    
  }
  
  ## --------------------------------------------------------
  
  ## set path back! 
  path = 'results/exome_genome'
  setwd(path)
  
}

## --------------------------------------------------------

sessionInfo()
