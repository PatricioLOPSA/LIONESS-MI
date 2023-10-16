# script para obtener los top 100 nodos de mayor fuerza en cada red
require(igraph)
require(dplyr)
require(data.table)
require(ggplot2)
require(stringr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lioness_output_dir <- commandArgs(trailingOnly = T)


  
  setwd(lioness_output_dir[1])
  
  file_ids <- list.files(pattern = "SSnet")
  nr_samples <- length(file_ids)
  
  output <- tibble(node = NA)
  
  #Lee SS_i, , la pasa a red de igraph, calcula strength
  for (i in 1:nr_samples) {
    
    g.net <- fread(file_ids[i], sep = "\t") %>% as.data.frame() %>% graph_from_data_frame(., directed=F)

    V(g.net)$strength <- strength(g.net, weights = E(g.net)$value)
    tcga_id <- str_remove_all(string = file_ids[i], pattern = "_a_min_q-complete.tsv-SSnet.tsv")
    
    node_centrality <- tibble(node = V(g.net)$name,
                                  strength = V(g.net)$strength) %>% as.data.frame()
    
    colnames(node_centrality) <- c("node", tcga_id)
    
    output <- merge(output, node_centrality, by ="node" ,all=T)
    output <- output[!is.na(output$node),]
    
    
    gc()
    
  }
  
  
output_filename <- paste0(lioness_output_dir[1],"_gene_targeting_score.tsv")  
  
fwrite(output, output_filename, sep = "\t")  
  

  
  



