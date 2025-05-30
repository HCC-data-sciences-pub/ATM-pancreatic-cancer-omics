
#########################################################
# Riyue Bao 06/17/2015
# This script combines expression matrix and tumor groups
# and run DEG detection.
# Update 08/01/2015
# Use raw counts as the input for limma-voom
#########################################################

rm(list=ls())

library(ggplot2)
library(scales)
library(ctc)
library(ape)
library(RColorBrewer)
library(reshape2)
library(limma)
library(edgeR)
library(Biobase)
library(VennDiagram)
library(ComplexHeatmap)
library(devtools)
# install_github('sinhrks/ggfortify')
library(ggfortify)

## http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
library(ggrepel)
library(circlize)
library(splines)

## https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
## https://github.com/cran/ggsci/blob/master/R/discrete-nejm.R
library(ggsci)
library(ggpubr)
library(plotrix)

#########################################################

run = 1 ## 1 - ABC vs D


# sample.ex = c()
sample.ex = c('T7290343', 'T5436736_Repeat')
sample.in = c()

## ------------------------------------------------------

cancer = "ATM_PAAD" 
fcs = c(2.0,1.5)
fdrs = c(0.05,0.10)
rawps = c()
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
                'mut' = '#000000', 'nonmut' = '#FFFFFF')

## ------------------------------------------------------

## paths
path = paste0('results/rnaseq/gene_expr')
setwd(path)

## ------------------------------------------------------

## files
group.file = '3.csv'
expr.file =  'ATM_PAAD.rnaseq.kallisto.raw.txi.coding.txt'
anno.file = 'gencode.v32.primary_assembly.annotation.maskPAR.gtf.geneinfo'

## ------------------------------------------------------------------------
## import data 
## ------------------------------------------------------------------------

print("Import data ...")

expr = read.delim(expr.file, header = T, stringsAsFactors = F, row.names = 1)
group = read.csv(group.file, header = T, stringsAsFactors = F, 
                   na.strings = c('na','NA'))
anno = read.delim(anno.file, header = T, stringsAsFactors = F)

# save(expr, file = paste0(expr.file,'.expr.RData'))
# load(paste0(expr.file,'.expr.RData'))

group$Group = paste0('C',group$cluster.x)

## ---------------------------------------------

dim(expr) ## 19784    40
dim(group) ## ## 40 107
dim(anno) ## 58395     6

table(group$Group)

group$Sample = group$sample

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

## ------------------------------------------------------------------------
## detect DEGs
## ------------------------------------------------------------------------

print("Detect DEGs ...")

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

if(run == 1) {
  
  group.col = 'Group'
  out.dir = file.path(out.dir,'Group')
  
} 


## ---------------------------------------------

## prepare design matrix and contrasts
table(group_sub[,group.col])
# mut nonmut 
# 6     16

Group = as.factor(group_sub[,group.col])

min.sm = NA
min.sm = min(table(Group))
if(!is.na(min.sm) & min.sm < 6) { min.sm = 6 }
# if(!is.na(min.sm) & min.sm < 6) { min.sm = 12 }
# if(!is.na(min.sm) & min.sm < 12) { min.sm = 12 } ## 02/18/2019: to set min.sm consistent through run=1/2/3....
print(paste0('Minimum sample per group = ', min.sm))

design = model.matrix(~ 0 + Group)

colnames(design) = gsub('Group','',gsub('Donor','',gsub('RunDate','',colnames(design))))
head(design)

contrast_matrix = makeContrasts(C1 - C2,
                                C1 - C3,
                                C1 - C4,
                                C2 - C3,
                                C2 - C4,
                                C3 - C4,
                                levels=design)

if(! dir.exists(out.dir)) { dir.create(out.dir, recursive = T) }

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
      
      ## from lei! Instead of using hard cpm cutoff...
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
            
            sample.anno = group_sub[,c('Sample',group.col)]
            sample.anno[is.na(sample.anno)] = 'NA'
            
            ## sort sample anno same as expression matrix
            sample.anno = sample.anno[order(match(sample.anno$Sample,
                                                  colnames(data.plot))),]
            row.names(sample.anno) = sample.anno$Sample
            # sample.anno = sample.anno[,-1]
            
            if(run == 1) {
              sample.anno.colors = list(
                Group = plot.colors
              )
            }
            
            plot.anno = HeatmapAnnotation(df = sample.anno[,! colnames(sample.anno)%in%
                                                             c('Sample')],
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
            
            if(run == 1) {
              sample.anno.colors = list(
                Group = plot.colors
              )
            }
            
            plot.anno = HeatmapAnnotation(df = sample.anno[,! colnames(sample.anno)%in%
                                                             c('Sample')],
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
    
    
    if(run == 1) {
      p5 = autoplot(data.pca, data = data.plot, colour = 'Group',
                    label = TRUE, label.size = 3)
    }
    
    pdf(file = paste0(output, '.',ncol(data.plot),'genes_',nrow(data.plot),'sm',
                      '.pca.pdf'), 
        width = 6.5, height = 5)
    print(p5)
    dev.off()
    
  }
}

## ----------------------------------------------------

## 01/07/2019
## check how many DEGs ovverlap with all samples or after excluding the 4 liver cell high UVM samples...

if(FALSE) {
  x = read.delim('C:/Users/rbaoc/Desktop/current_proj/lukeProj/UvealMelanoma_grant/results/rnaseq/diff_expr/Group/limmavoomweighted/JLUVM.19sm.coding.cpm3.limmavoomweighted.YvsN.txt',stringsAsFactors = F)
  y = read.delim('C:/Users/rbaoc/Desktop/current_proj/lukeProj/UvealMelanoma_grant/results/rnaseq/diff_expr/Group/limmavoomweighted/JLUVM.15sm.coding.cpm3.limmavoomweighted.YvsN.txt',stringsAsFactors = F)
  x = x[x$P.Value<0.005 & abs(x$FC)>=2,]
  dim(x) ## 76 genes from 19sm analysis
  x = merge(x, y, by = 'Gene', all.x = T)
  dim(x) ## 16
  write.csv(x, file = paste0(out.dir,'/limmavoomweighted/19sm.vs.15sm.DEGs.csv'),row.names = F)

  
  y = y[y$P.Value<0.005 & abs(y$FC)>=2,] 
  dim(y) ## 60
  
  length(intersect(x$Gene,y$Gene)) ## 42
  
  42/76 ## 0.5526316
  42/60 ## 0.7
}

#########################################################

print("...Done!")

sessionInfo()


