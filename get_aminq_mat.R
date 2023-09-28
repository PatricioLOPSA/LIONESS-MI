library(data.table)
library(dplyr)

sysargs <- commandArgs(trailingOnly = T)

#expression matrix: all samples
exp <- fread(sysargs) %>% as.data.frame()



ids <- colnames(exp[,-1])
ids <- paste(ids, "_a_min_q.tsv", sep = "")


for (i in 1:length(ids)) {

  ss_exp <- exp[,-(i+1)]

  fwrite(ss_exp, file = ids[i], sep = "\t", )
}


