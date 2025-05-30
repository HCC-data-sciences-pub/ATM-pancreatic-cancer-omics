

rm(list = ls())

## -------------------------------------------

# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")

library(immunedeconv)
library(PerformanceAnalytics)
library(ComplexHeatmap)
library(circlize)
library(broom)

## -------------------------------------------

## settings
if(TRUE) {
  methods = c('quantiseq', 'timer', 'cibersort', 'cibersort_abs', 'mcp_counter', 'xcell', 'epic')
  my.cancer = 'paad'
  
  plot.colors = c('atm_paad' = '#CCCC00',
                  'tcga_paad' = '#00CCCC')
}

## -------------------------------------------

path = 'results/rnaseq/gene_expr'
setwd(path)

## -------------------------------------------

## import data and select samples 
if(TRUE) {
  expr.file = 'ATM_PAAD_tpm.TCGA_PAAD_rsem.rds'
  
  expr = readRDS(expr.file)
  
}

## -------------------------------------------

## prep inputs
if(TRUE) {
  
  expr = as.matrix(expr)
  row.names(expr) = gsub('[!]\\S+$','', row.names(expr))
  dim(expr) ## 18427   205
  
  write.csv(expr,
            file = paste0(expr.file,
                          '.sm',ncol(expr),'.immunedeconv.input.csv'))
  
  ## as input for cibersort online ... it does NOT take csv files
  write.table(expr,
              file = paste0(expr.file,
                            '.sm',ncol(expr),'.immunedeconv.input.txt'),
              sep = '\t',col.names = NA,row.names = T, quote = F)
}

## -------------------------------------------

## run deconv
if(TRUE) {
  deconv.out = list(
    quantiseq = NULL,
    timer = NULL ,
    cibersort = NULL,
    cibersort_abs = NULL ,
    mcp_counter = NULL ,
    xcell = NULL ,
    epic = NULL 
  )
  
  deconv.out.df = NULL 
  
  for(my.m in methods) {
    print(my.m)
    my.df = NULL 
    
    if(my.m %in% c('cibersort','cibersort_abs')) {
      
      # next 
      set_cibersort_binary("data/cibersort/CIBERSORT.R")
      set_cibersort_mat("data/cibersort/LM22.txt")
      
      my.df = deconvolute(expr, my.m)
      
    } else if (my.m == 'timer') {
      my.df = deconvolute(expr, "timer",
                          indications=rep(my.cancer,ncol(expr)))
    } else {
      my.df = deconvolute(expr, my.m)
    }
    
    deconv.out[[my.m]] = my.df
    
    deconv.out.df = rbind(deconv.out.df,
                         data.frame(
                           method = my.m,
                           my.df,
                           stringsAsFactors = F
                         ))
    
  }
  
  saveRDS(deconv.out,
          paste0(expr.file,
                 '.sm',ncol(expr),'.immunedeconv.out.rds'))
  
  write.csv(deconv.out.df,
            paste0(expr.file,
                   '.sm',ncol(expr),'.immunedeconv.out.csv'))
}

## -------------------------------------------

## correlation between different methods for CD8
if(TRUE) {
  
  ## read in files 
  if(TRUE) {
    deconv.out = readRDS(paste0(expr.file,
                                '.sm',ncol(expr),'.immunedeconv.out.rds'))
    deconv.out.df = read.csv(paste0(expr.file,
                                    '.sm',ncol(expr),'.immunedeconv.out.csv'),
                             stringsAsFactors = F, row.names = 1)
    
    samples = deconv.out.df[-c(1:2)]
    
    deconv.out.df[grep('CD8', deconv.out.df$cell_type),1:2,drop=F]
    # method                   cell_type
    # 8    quantiseq                 T cell CD8+
    # 14       timer                 T cell CD8+
    # 19 mcp_counter                 T cell CD8+
    # 36       xcell           T cell CD8+ naive
    # 37       xcell                 T cell CD8+
    # 38       xcell  T cell CD8+ central memory
    # 39       xcell T cell CD8+ effector memory
    # 71        epic                 T cell CD8+
    
    x = data.frame(table(deconv.out.df$cell_type))
    x[grep('T cell',x$Var1),]
    x[grep('NK',x$Var1),]
  }

  ## -------------------------------------------
  
  if(TRUE) {
    
    my.cell = 'T cell CD8+'
    # my.cell = 'NK cell'

    rm(data.plot)
    data.plot = deconv.out.df[deconv.out.df$cell_type==my.cell,,drop=F]
    
    dim(data.plot)
    
    row.names(data.plot) = data.plot$method
    data.plot = data.plot[,colnames(samples)]
    dim(data.plot) ## 5 22
    
    ## add cibersort results ...
    ## I manually ran this on cibersort website ...
    ## input file = ECOG.rnaseq_encode28_coding_tpm.csv.sm22.immunedeconv.input.txt
    if(FALSE) {
      
      cibersort = read.csv('ECOG.rnaseq_encode28_coding_tpm.csv.sm22.immunedeconv.input.txt.cibersort_rel.csv',stringsAsFactors = F, row.names = 1)
      dim(cibersort)
      cibersort = data.frame(t(cibersort))
      
      cibersort = cibersort[,colnames(samples)]
      
      cibersort_abs = read.csv('ECOG.rnaseq_encode28_coding_tpm.csv.sm22.immunedeconv.input.txt.cibersort_abs.csv',stringsAsFactors = F, row.names = 1)
      dim(cibersort_abs)
      cibersort_abs = data.frame(t(cibersort_abs))
      
      cibersort = cibersort[,colnames(samples)]
      
      rm(row.select)
      if(my.cell == 'my.cell') {
        
        row.select = row.names(cibersort)[grep('CD8',row.names(cibersort))]
        ## "T.cells.CD8"
      } else if (my.cell == 'NK cell') {
        
        row.select = row.names(cibersort)[grep('NK',row.names(cibersort))]
        ## "NK.cells.resting"   "NK.cells.activated"
      }
      
      print(row.select)
      
      rm(x,y)
      x = cibersort[row.names(cibersort) %in% row.select,,drop=F]
      y = cibersort_abs[row.names(cibersort_abs) %in% row.select,,drop=F]
      
      if(length(row.select)==1) {
        row.names(x) = 'cibersort'
        row.names(y) = 'cibersort_abs'
      } else {
        row.names(x) = paste0('cibersort',1:length(row.select))
        row.names(y) = paste0('cibersort_abs',1:length(row.select))
      }
      
      data.plot = rbind(data.plot,
                        x,
                        y)
      
    }
    
    data.plot[,1:5]
    data.plot = t(data.plot)
    
    write.csv(data.plot,
              paste0(expr.file,
                     '.sm',ncol(expr),'.immunedeconv.out.',my.cell,'.cor.csv'))
    
    ## -------------------------------------------
    
    pdf(file = paste0(expr.file,
                      '.sm',ncol(expr),'.immunedeconv.out.',my.cell,
                      '.cor.pearson.pdf'),
        width = 6, height = 6)
    chart.Correlation(data.plot, histogram=TRUE, pch=19,
                      method = 'pearson')
    dev.off()
    
    pdf(file = paste0(expr.file,
                      '.sm',ncol(expr),'.immunedeconv.out.',my.cell,
                      '.cor.spearman.pdf'),
        width = 6, height = 6)
    chart.Correlation(data.plot, histogram=TRUE, pch=19,
                      method = 'spearman')
    dev.off()
    
    ## -------------------------------------------
    
    pdf(file = paste0(expr.file,
                      '.sm',ncol(expr),'.immunedeconv.out.',my.cell,
                      '.cor.pearson.exEPIC.pdf'),
        width = 6, height = 6)
    chart.Correlation(data.plot[,!colnames(data.plot) %in% c('epic')], 
                      histogram=TRUE, pch=19,
                      method = 'pearson')
    dev.off()
    
    pdf(file = paste0(expr.file,
                      '.sm',ncol(expr),'.immunedeconv.out.',my.cell,
                      '.cor.spearman.exEPIC.pdf'),
        width = 6, height = 6)
    chart.Correlation(data.plot[,!colnames(data.plot) %in% c('epic')], 
                      histogram=TRUE, pch=19,
                      method = 'spearman')
    dev.off()
    
  }

}

## -------------------------------------------

## heatmap 
if(TRUE) {
  
  deconv.out = readRDS(paste0(expr.file,
                             '.sm',ncol(expr),'.immunedeconv.out.rds'))
  deconv.out.df = read.csv(paste0(expr.file,
                                  '.sm',ncol(expr),'.immunedeconv.out.csv'),
                           stringsAsFactors = F, row.names = 1)
  
  my.m = 'xcell' ## xcell cibersort
  
  ## atm paad samples 
  if(TRUE) {
    
    ## prep inputs 
    if(TRUE) {
      data.anno = NULL 
      data.anno = data.frame(Sample = colnames(deconv.out.df)[-c(1:2)],
                             stringsAsFactors = F)
      data.anno$Cohort = NA 
      data.anno$Cohort[grep('^T\\d+',data.anno$Sample)] = 'atm_paad'
      data.anno$Cohort[grep('^TCGA',data.anno$Sample)] = 'tcga_paad'
      
      table(data.anno$Cohort)
      # atm_paad tcga_paad 
      # 22       183
      
      ## -------------------------------------------
      
      data.plot = NULL 
      data.plot = deconv.out.df[deconv.out.df$method==my.m,]
      row.names(data.plot) = data.plot$cell_type
      data.plot =  data.plot[,data.anno$Sample]
      # data.plot =  data.plot[,data.anno$Sample[data.anno$Cohort=='atm_paad']]
      
      dim(data.plot)
      
    }

    ## -------------------------------------------

    ## plot heatmap
    if(TRUE) {

      centered_data = as.matrix(data.plot)
      
      if(TRUE) {
        
        ## add annotation 
        
        sample.anno = data.anno
        sample.anno[is.na(sample.anno)] = 'NA'

        ## sort sample anno same as expression matrix
        sample.anno = sample.anno[sample.anno$Sample %in% colnames(data.plot),]
        sample.anno = sample.anno[order(match(sample.anno$Sample,
                                              colnames(data.plot))),]
        row.names(sample.anno) = sample.anno$Sample
        # sample.anno = sample.anno[,-1]
        print(all.equal(sample.anno$Sample, colnames(data.plot)))
        
        sample.anno.colors = list(
          Cohort = plot.colors
        )
        
        plot.anno = HeatmapAnnotation(df = sample.anno[,c('Cohort'),drop=F],
                                      col = sample.anno.colors)
        
      }
      
      # myheatcol = colorRamp2(c(-1.6, 0, 1.6), c("#6060a8", "#FFFFFF", "#f01830"))
      if(my.m %in% c('xcell','cibersort')) {
        myheatcol = colorRamp2(c(0, 0.5), c("white", "red"))
      }
      
      
      col.title = paste0(expr.file, ' ',ncol(centered_data), 
                         ' samples')
      row.title = paste0(nrow(centered_data), ' genes,')
      
      p1 = Heatmap(centered_data,
                   na_col = "#000000",
                   col = myheatcol,
                   rect_gp = gpar(col = NA),
                   show_heatmap_legend = T,
                   column_title = col.title,
                   row_title = row.title,
                   # column_title_side = 'bottom',
                   column_names_side = 'bottom',
                   row_dend_width = unit(5, "cm"),
                   column_dend_height = unit(5, "cm"),
                   # km = 2,
                   cluster_rows = T,
                   cluster_columns = T,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "ward.D2",
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "ward.D2",
                   show_row_names = T,
                   top_annotation = plot.anno,
                   heatmap_legend_param = list(title = 'log2(Expression)', 
                                               color_bar = "continuous")
      )
      
      p2 = Heatmap(centered_data,
                   na_col = "#000000",
                   col = myheatcol,
                   rect_gp = gpar(col = NA),
                   show_heatmap_legend = T,
                   column_title = col.title,
                   row_title = row.title,
                   # column_title_side = 'bottom',
                   column_names_side = 'bottom',
                   row_dend_width = unit(5, "cm"),
                   column_dend_height = unit(5, "cm"),
                   # km = 2,
                   cluster_rows = T,
                   cluster_columns = F,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "ward.D2",
                   # clustering_distance_columns = "euclidean",
                   # clustering_method_columns = "ward.D2",
                   show_row_names = T,
                   top_annotation = plot.anno,
                   heatmap_legend_param = list(title = 'log2(Expression)', 
                                               color_bar = "continuous")
      )

      output = paste0(expr.file,'.',ncol(centered_data),'sm_',nrow(centered_data),
                      'genes.',my.m,'.heatmap')
      pdf(paste0(output,'.pdf'), 
          width = 16, height = 12)
      print(p1)
      dev.off()
      
      pdf(paste0(output,'.srtBySample.pdf'), 
          width = 16, height = 10)
      print(p2)
      dev.off()
      
      write.table(data.plot,
                  file = paste0(output,'.input'),
                  col.names = NA, row.names = T, sep = '\t', quote = F)
      
      ## -------------------------------------------
      
      # ## get row order...
      # data.plot = deg.sig[order(match(gsub('[!]\\w+$','',deg.sig$Gene),
      #                                 row.names(centered_data)[row_order(p8)])),]
      # 
      # write.table(data.plot,
      #             file = paste0(output,'.row_order'),
      #             col.names = T, row.names = F, sep = '\t', quote = F)
      # 
      # ## print p-value heatmap...
      # data.plot = data.plot[,c('Gene','adj.P.Val')]
      # row.names(data.plot) = gsub('[!]\\w+$','',data.plot$Gene)
      # data.plot = data.plot[,-1]
      # data.plot[data.plot>=0.10] = 1
      # data.plot[data.plot<0.10 & data.plot>=0.05] = 0.10
      # data.plot[data.plot<0.05 & data.plot>=0.01] = 0.05
      # data.plot[data.plot<0.01] = 0.01
      # data.plot = as.matrix(data.plot)
      # 
      # myheatcol = colorRamp2(c(0, 0.01,0.05,0.10,1), 
      #                        c("#000000", "#606060", "#A0A0A0","#E0E0E0","#FFFFFF"))
      # # myheatcol = colorRamp2(c(-1.6, 0, 1.6), c("#6060a8", "#FFFFFF", "#f01830"))
      # 
      # col.title = paste0(expr.file, ' ',ncol(centered_data), 
      #                    ' samples')
      # row.title = paste0(nrow(centered_data), ' genes, FDR<',fdr, ' fc',fc)
      # 
      # p10 = Heatmap(data.plot,
      #               na_col = "#000000",
      #               col = myheatcol,
      #               rect_gp = gpar(col = NA),
      #               show_heatmap_legend = T,
      #               column_title = col.title,
      #               row_title = row.title,
      #               # column_title_side = 'bottom',
      #               column_names_side = 'bottom',
      #               row_dend_width = unit(5, "cm"),
      #               column_dend_height = unit(5, "cm"),
      #               # km = 2,
      #               cluster_rows = F,
      #               cluster_columns = F,
      #               clustering_distance_rows = "euclidean",
      #               clustering_method_rows = "ward.D2",
      #               clustering_distance_columns = "euclidean",
      #               clustering_method_columns = "ward.D2",
      #               show_row_names = T,
      #               # top_annotation = plot.anno,
      #               heatmap_legend_param = list(title = 'p-value')
      # )
      # 
      # pdf(paste0(output,'.row_order.pvalue.pdf'), 
      #     width = 6, height = 8)
      # print(p10)
      # dev.off()
      
    }
    
    ## -------------------------------------------
    
    ## compare atm_paad vs tcga_paad
    if(TRUE) {
      
      dim(data.plot)
      
      data.stats = NULL       
      for(i in 1:nrow(data.plot)) {
        my.gene = row.names(data.plot)[i]
        print(my.gene)
        
        my.df = NULL 
        my.df = data.frame(t(data.plot[i,,drop=F]),
                           stringsAsFactors = F)
        my.df = data.frame(Sample = row.names(my.df),
                           my.df,
                           stringsAsFactors = F)
        
        my.df = merge(my.df, data.anno, by = 'Sample')
        rm(x,y, z)
        x = my.df[my.df$Cohort=='atm_paad',2]
        y = my.df[my.df$Cohort=='tcga_paad',2]
        z = tidy(wilcox.test(x,y))
        
        data.stats  = rbind(data.stats,
                            data.frame(Gene = my.gene,
                                       mean.1 = mean(x),
                                       mean.2 = mean(y),
                                       z,
                                       stringsAsFactors = F))
        
      }
      
      data.stats$p.adj = p.adjust(data.stats$p.value, method = 'fdr')
      
      write.csv(data.stats,
                file = paste0(expr.file,'.',ncol(centered_data),'sm_',nrow(centered_data),
                              'genes.',my.m,'.stats.csv'))
      
    }
    
  }
  
}

## -------------------------------------------

sessionInfo()
