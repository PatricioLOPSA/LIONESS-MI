library(dplyr)
library(reshape2)
library(data.table)

#arg 1 = lioness output nets file
#arg 2 = new output dir name
#arg 3 = number of edges

sysargs <- commandArgs(trailingOnly = T)

dir_nets <- sysargs[1]
dir_output <- sysargs[2]
edges_number <- as.numeric(sysargs[3])

launch <- getwd()

abs_dir_nets <- paste(launch,dir_nets,sep = "/")
abs_dir_output <- paste(launch, dir_output, sep = "/")

dir.create(dir_output)


setwd(dir_nets)




#Get SS net file ids------------------------------------------------------------
#Leer ids de cada red en el dir-------------------------------------------------
file_ids <- list.files() 
file_ids <- file_ids[ !(file_ids =="alpha_net.tsv")]
nr_samples <- length(file_ids) 
setwd("..")

for (i in 1:nr_samples) {
  
  setwd(abs_dir_nets)
  SS_net <- fread(file_ids[i], sep = "\t") %>% as.data.frame()
  #SS_net <- SS_net[SS_net$value > 0,]
  SS_net <- head(SS_net,edges_number)
  setwd(abs_dir_output)
  fwrite(SS_net, file_ids[i], sep = "\t")
  
}

#move alphanet
alpha_abs_dir <- paste(abs_dir_nets,"alpha_net.tsv", sep = "/")

file.copy(from = alpha_abs_dir,
          to = abs_dir_output)




# #Quitar partes inecesarias del id resultantes del pipeline----------------------
# tcga_ids <- str_remove_all(string = file_ids, 
#                            pattern = "_a_min_q-complete.tsv-100k-SSnet.tsv")
# 
