

## prep gistic2 input ...

library(tidyr)

datatype = 'genome' ## exome or genome 
project = 'ATM_PAAD'
aligner = 'bwamem'
caller = 'freec'
out = paste0(project, '.',datatype, '.', aligner, '.',caller)

path = paste0('/results/',datatype)

setwd(path)

in.dir = paste0(caller)
ratio.files = list.files(path = in.dir, pattern = "*.bam_ratio.txt$",
                       recursive = T)

ratio.files
length(ratio.files)

ratio.combined = NULL 

## gistic2 values = log2(depth ratio *2) - 1
for(my.file in ratio.files) {
  
  print(my.file)
  my.df = NULL 
  my.sample = NULL 
  
  my.sample = gsub('[/]\\S+','',my.file)
  
  my.df = read.delim(file.path(in.dir, my.file),
                     stringsAsFactors = F)
  print(dim(my.df))
  
  my.df = my.df[!is.na(my.df$MedianRatio) & my.df$MedianRatio!=-1,]
  print(dim(my.df))
  
  my.df$gistic2.input = log2(my.df$MedianRatio * 2)-1
  
  my.df = my.df %>% separate(Gene, sep = '[:-]', 
                             c('chr','start','end'))
  
  ratio.combined = rbind(ratio.combined,
                         data.frame(Sample = my.sample,
                         my.df[,c('chr','start','end')],
                         markers = 100,
                         gistic2.input = my.df$gistic2.input,
                         stringsAsFactors = F)
  )
  
}

head(ratio.combined)

save(ratio.combined, 
     file = paste0(out,'.gistic2.seg.RData'))

write.table(ratio.combined,
            file = paste0(out,'.gistic2.seg.txt'),
            sep = '\t', col.names = F, row.names = F, quote = F)