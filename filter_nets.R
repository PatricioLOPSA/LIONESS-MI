library(data.table, quietly = F)
library(dplyr, quietly = T)


 sysargs <- commandArgs(trailingOnly = T)
 
 setwd(sysargs[1])

#setwd("/datos/ot/plos/aracne-multicore/launch/Luadtp-SSnets-stgI")

files <- list.files(pattern="SSnet.tsv") %>% as.list()

#aquí empieza la function:
filter_nets <- function(x){
net_id <- as.character(x)
net <- fread(file = net_id) %>% as.data.frame()
net <- net[order(- net[,3]),]
nr_edges <- as.numeric(sysargs[2])
#nr_edges <- 10000

SS_name <- paste("filtered-",net_id, sep = "")

fwrite(head(net,nr_edges), SS_name,sep = "\t")
}

lapply(files, filter_nets)


