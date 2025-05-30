
rm(list = ls())

library(ConsensusClusterPlus)

# rm(list=ls(all=TRUE))
# ls()

# args <- commandArgs(TRUE)
# file = args[1]
# out_dir = args[2]
# type = args[3]
# print(file)
# print(out_dir)
# print(type)

# rnaseq.ex = c('T7290343')
rnaseq.ex = c('T7290343', '5436736_Repeat')

file = "ATM_PAAD.rnaseq.kallisto.raw.txi.coding.exT7290343_T5436736_Repeat.TMM.logCPM.sm22.csv"
out_dir = "ccp_output"

# file = "NBL.sm161.rnaseq.kallisto.raw.txi.coding.voom.tcellgenes.matrix"
# out_dir = "ccp_output_tcellgenes"


path="results/rnaseq/gene_expr"
setwd(path)

## global parameters"
type = "sm"
# type = "gene"
max_k = 20 #20
reps = 2000 #2000
pct_item = 0.8 ## sample
pct_feature = 1 ## gene
clr_method = "hc"
clr_dist_gene = "euclidean"
clr_linkage_gene = "ward.2D"
# clr_dist_sm = "pearson"
# clr_linkage_sm = "average"
clr_dist_sm = "euclidean"
clr_linkage_sm = "ward.D2" 
seed = 119
plot_type = "pdf"
colors=rev(heat.colors(10))

# setwd(path)

## prepare input data

data = read.table(file, header=T,sep="\t",na.strings = "NA",row.names=1)
#data = data[,1:334] ## after re-importation, remove the last col "mad"

## calculate MAD for each row, and manually filter the resulting matrix
# mads = apply(data,1,mad,na.rm=T)
# head(mads)
# length(mads)
# str(mads)
# data_mads = cind(data,mads)
# write.table(cbind(data,mads),paste0(file,".srtByMad"),sep="\t",col.names=T,row.names=T,quote=F)

## normalize by row using Pearson correlation distance
## sample cluster
if(type=="sm")
{
  ## listwise deletion of missing; if turned on, gene rows may reduce a lot, which led to wierd CCP results
  ## In general, cr is better for SM clustering. However, if your data has lots of NAs (e.g. GDAC rnaseq), better use eu distance 
  if(clr_dist_sm=="pearson") { data = na.omit(data) } 
  
  centered_data = as.matrix(sweep(data,1,apply(data,1,median,na.rm=T)))
  results = ConsensusClusterPlus(centered_data,maxK=max_k,reps=reps,pItem=pct_item,pFeature=pct_feature,title=out_dir,clusterAlg=clr_method,distance=clr_dist_sm,innerLinkage=clr_linkage_sm,finalLinkage=clr_linkage_sm,seed=seed,plot=plot_type,tmyPal=colors,writeTable=T,verbose=T)
  save(results,file=paste0(file,".ccp_results.sm.RData"))
}
## gene cluster
if (type=="gene") 
{ 
  centered_data = t(scale(t(data), scale=F))
  #centered_data = scale(data, scale=F) ## scale by column for row clustering - the clustering result is horrible, stick to the original plan
  centered_data_t = t(centered_data) 
  results = ConsensusClusterPlus(centered_data_t,maxK=max_k,reps=reps,pItem=pct_item,pFeature=pct_feature,title=out_dir,clusterAlg=clr_method,distance=clr_dist_gene,innerLinkage=clr_linkage_gene,finalLinkage=clr_linkage_gene,seed=seed,plot=plot_type,tmyPal=colors,writeTable=T,verbose=T)
  save(results,file=paste0(file,".ccp_results.RData"))
} 

# class(results)
# str(results)

