require(CoDiNA)
require(data.table)
require(dplyr)
require(stringr)

# dir_exp <- "/datos/ot/plos/aracne-multicore/launch/" ###***
# 
# exp_file <- "mock_dataset.tsv"
# dir_nets <- "mock"

sysargs <- commandArgs(trailingOnly = T)

dir_nets <- sysargs[1]
exp_file <- sysargs[2]



#Conseguir todos los genes dentro de la matriz de expresión--------------------
#setwd(dir_exp)
exp_mat <- fread(exp_file) %>% as.data.frame()
node_list <- exp_mat[,1]

sample_id <- colnames(exp_mat)
faulty_id <- sample_id[2] #quitar cuando se corrija la duplicación de red 1



empty_nodes <- data_frame(Var1=node_list,Var2=node_list,value=0)





#Get SS net file ids------------------------------------------------------------
setwd(dir_nets)
file_ids <- list.files() 

faulty_file <- list.files(, pattern = faulty_id) #quitar cuando se corriga la duplicación de red 1
file_ids<- file_ids[!(file_ids %in% faulty_file)] #quitar cuando se corriga la duplicación de red 1

file_ids <- file_ids[!(file_ids %in% "alpha_net.tsv")]


#Format sample Ids -------------------------------------------------------------
nr_samples <- length(file_ids) 

tcga_ids <- str_remove_all(string = file_ids, 
                           pattern = "_a_min_q-complete.tsv-100k-SSnet.tsv")
tcga_ids_mods <- gsub("-","_",tcga_ids)


# index <- as.character(1:nr_samples)
# n_index <- paste("n",index,sep = "")
# net_ids <- data_frame(tcga_id=tcga_ids, net_index=n_index)



#Create CoDiNA input list of networks-------------------------------------------
input <- c()

for (i in 1:nr_samples) {
  
  net <- fread(file_ids[i], sep = "\t") %>% as.data.frame()
  
  #La siguiente linea es en caso de que se quiera omitir el filtrado de CoDiNA
  net$value <- ifelse(net$value >0, 1,-1)
  
  net <- rbind(net, empty_nodes)
  input <- c(input, list(net))
  
}


#Run CoDiNA analsis-----------------------------------------------------------------
diffnet <- MakeDiffNet(input, Code = tcga_ids_mods) %>% as.data.frame()

# Filter out CoDiNA stats column except phi tilde and get edge frequencies------
phi_tilde <- diffnet$Phi_tilde
reg <- diffnet$Node.1
tar <- diffnet$Node.2

codina_stats <- c("Phi","Group","Score_center","Node.1","Node.2","Phi_tilde",
                  "Score_Phi","Score_Phi_tilde","Score_internal","Score_ratio")

#Filter out non edge-value columns
diffnet <- diffnet[,!(names(diffnet) %in% codina_stats)]

#Sum by row all edge values different from 0 to get edge frequencies
row_counts <- rowSums(diffnet != 0)

#bind reg, tar, edge values, edge frequencies and Phi tilde, and arrange by edge frequencies
diffnet <- cbind.data.frame(reg, tar, diffnet, row_counts,phi_tilde) %>% 
  arrange(desc(row_counts))

#save output in same directory as expression matrix-----------------------------
setwd("..")
output_name <- paste(dir_nets,"CoDiNA_diffnet.tsv",sep = "-")
fwrite(diffnet, output_name, sep = "\t")







