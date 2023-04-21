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
SS <- nr_samples * (alpha_net - a_min_q) + a_min_q

elify <- function(net){
	    require(reshape2)
    net[lower.tri(net, diag = T)] <- NA
        net <- melt(net)
        net <- net[!is.na(net$value),]
	    return(net)
	  }


SS <- elify(SS)

SS <- SS[order(- SS[,3]),]


SS_name <- paste(sysargs[2], "-100k-SSnet.tsv", sep = "")
fwrite(head(SS,100000), SS_name,sep = "\t")
