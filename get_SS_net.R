library(data.table, quietly = T)
library(dplyr, quietly = T)

sysargs <- commandArgs(trailingOnly = T)


agg_exp <- fread(file = sysargs[1])
nr_samples <- length(agg_exp) - 1

alpha_net <- fread("alpha_net.tsv") %>% as.matrix()
a_min_q <- fread(file = sysargs[2]) %>% as.matrix()

rownames(a_min_q) <- colnames(a_min_q)
rownames(alpha_net) <- colnames(alpha_net)

#Lioness equation
SS <-round( nr_samples * (alpha_net - a_min_q) + a_min_q,3)

elify <- function(net){
    require(reshape2)
    net[lower.tri(net, diag = T)] <- NA
    net <- melt(net)
    net <- net[!is.na(net$value),]
    return(net)
  }


SS <- elify(SS)


SS_name <- paste(sysargs[2], "-SSnet.tsv", sep = "")
fwrite(SS, SS_name,sep = "\t")




#table <- fread(sysargs)

#print(table[1:5,1:5])

#library(data.table)
# adj <- fread("luadtpIII_aminq_exp-complete.tsv")  %>% as.matrix() #, col_names = T, delim = "\t") %>% as.data.frame()
# adj2 <- adj
# adj3 <- adj + adj2
# adj3[1:5,1:5]
# toy_data <- exp[1:10,1:10]
# fwrite(toy_data, file = "toy_data.tsv", sep = "\t")
                                                              
