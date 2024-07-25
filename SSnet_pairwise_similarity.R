library(dplyr)
library(reshape2)
library(data.table)
library(SummarizedExperiment)
library(stringr)


#Arg 1 = SSNs directory
#Arg 2 = output tsv file name

sysargs <- commandArgs(trailingOnly = T)

#set wd
dir <- sysargs[1]
setwd(dir)

#Read all SSN ids-------------------------------------------------
file_ids <- list.files() 
file_ids <- file_ids[! (file_ids =="alpha_net.tsv")]

#remove process identifiers from the network building process---------------------------
tcga_ids <- str_remove_all(string = file_ids, pattern = "_a_min_q-complete.tsv-SSnet.tsv")
nr_samples <- length(file_ids) 

#Setup parallel backend---------------------------------------------------------
require(parallel)
require(foreach)
require(doParallel)

cores <-  40

my.cluster <- parallel::makeCluster(cores, type = "FORK")
registerDoParallel(cl = my.cluster)


SSN_pwise_sim <- foreach(i = 1:nr_samples, .combine = "cbind") %dopar% {
 SS_i <- fread(file_ids[i]) %>% as.data.frame()
 SS_i$aux <- paste(SS_i[,1], SS_i[,2], sep = "-")
 
 outputaux <- matrix(NA, nrow = nr_samples, ncol = 1)
 
 #Read all nets to compare against SS_i
 for (j in 1:nr_samples) {
   SS_j <- fread(file_ids[j]) %>% as.data.frame()
   SS_j$aux <- paste(SS_j[,1], SS_j[,2], sep = "-")
   
   #get intersect of edges between SS_i y SS_j and save in [i,j]
   
   outputaux[j] <- length(intersect(SS_i$aux, SS_j$aux))
   
 }

 return(outputaux)
  
}

stopCluster(cl = my.cluster)
#-------------------------------------------------------------------------------

colnames(SSN_pwise_sim) <- tcga_ids
rownames(SSN_pwise_sim) <- tcga_ids

setwd("..")
fwrite(SSN_pwise_sim, file = sysargs[2], sep = "\t")



