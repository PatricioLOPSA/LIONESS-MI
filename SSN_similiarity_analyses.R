library(dplyr)
library(CoDiNA)
library(igraph)
library(data.table)
library(stringr)
library(SummarizedExperiment)
library(xCell)
library(biomaRt)
library(janitor)
library(ggsurvfit)
library(survival)
library(survminer)
library(RColorBrewer)
library(pheatmap)
library(rstatix)
library(corrplot)
library(ggsankey)

#script for patient similarity graphs/matrix analyses:
#graph-clustering with louvain, CoDiNA analysis for getting edge frequencies within clusters,
#Xcell Cell-type enrichment and Clinical Overrepresentation testing

#Define directories

#Location of FILTERED SSNs
ssn_dir <- "/datos/ot/plos/aracne-multicore/launch/filtered_10k_SSnets_surv_luadtp"

#Location of previously obtained Similarity Matrix
dir <- "/datos/ot/plos/aracne-multicore/launch"

#Location of clinical data object
clinic_dir <- "/datos/ot/plos/"

#Location of dir to store Codina tables
codina_outputs <- "codina_outputs"



#Define  functions:
elify <- function(net){
  require(reshape2)
  net[lower.tri(net, diag = T)] <- NA
  net <- melt(net)
  net <- net[!is.na(net$value),]
  return(net)
}
cooler_row_merge <- function(x,y){
  
  a <- merge(x,y, by=0)
  rownames(a) <- a$Row.names
  a$Row.names <- NULL
  return(a)
  
  
}
matches_pattern <- function(element, patterns) {
  any(sapply(patterns, function(pattern) grepl(pattern, element)))
}


#1.- import similarity matrices, expression matrix and clinical data.-----------

setwd(dir)

#expression and shared number of edges per network matrix
expr_matrix <- fread("survival_luadtp_matrix.tsv") %>% as.data.frame() 
sim_matrix <- fread("filtered10k_surv_SSnets_similarity.tsv", sep = "\t") %>% as.data.frame()
title <- "Pairwise single sample network similarity -10k networks"


sample_id <- colnames(expr_matrix)
faulty <- sample_id[2]

rownames(sim_matrix) <- colnames(sim_matrix)



sim_matrix <- sim_matrix[! rownames(sim_matrix) == faulty, ! colnames(sim_matrix) == faulty]

file_ids <- colnames(sim_matrix)
tcga_ids <- str_remove_all(string = file_ids, 
                           pattern = "_a_min_q-complete.tsv-SSnet.tsv")
tcga_ids <- str_remove_all(string = tcga_ids, 
                           pattern = "filtered-")
colnames(sim_matrix) <- tcga_ids
rownames(sim_matrix) <- tcga_ids

#transform similarity matrix form # of edges to jaccard similarity index
n_edges <- 10000
jaccard_imat <- sim_matrix/((n_edges*2)-sim_matrix)
sim_matrix <- jaccard_imat
rm(jaccard_imat)

#import clinical data:
setwd(clinic_dir)
SE <- "rnaseq-luad-clinicstg.rds"
clinic_info <- readRDS(SE) 
extract_clinic <- function(exp_mat){
  clinic_df <- colData(clinic_info) %>% as.data.frame() %>% filter(., rownames(.) %in% names(exp_mat))
  return(clinic_df)
}
clinic <- extract_clinic(expr_matrix)

#Format for enrichment and clustering analyses
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage I"] <- "stage1"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IA"] <- "stage1"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IB"] <- "stage1"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage II"] <- "stage2"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IIA"] <- "stage2"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IIB"] <- "stage2"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage III"] <- "stage3"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IIIA"] <- "stage3"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IIIB"] <- "stage3"
clinic$ajcc_pathologic_stage[clinic$ajcc_pathologic_stage=="Stage IV"] <- "stage4"
clinic$ajcc_pathologic_stage <- as.factor(clinic$ajcc_pathologic_stage)


clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T1"] <- "T1"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T1a"] <- "T1"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T1b"] <- "T1"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T2"] <- "T2"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T2a"] <- "T2"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T2b"] <- "T2"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T3"] <- "T3"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="T4"] <- "T4"
clinic$ajcc_pathologic_t[clinic$ajcc_pathologic_t=="TX"] <- "TX"
clinic$ajcc_pathologic_t <- as.factor(clinic$ajcc_pathologic_t)





#flatten similarity matrix to create patient-patient edge list
flat_matrix <- elify(as.matrix(sim_matrix))
flat_matrix <- flat_matrix[order(- flat_matrix[,3]),]
head(flat_matrix)
#viz jaccard indices distribution
plot(density(flat_matrix[,3]), ylim=c(0,1))
#-------------------------------------------------------------------------------


#Cluster patient similiarity matrices:------------------------------------------
#Generate clusters with graph based methods: Louvain
patient_graph <- graph_from_data_frame(flat_matrix, directed = F)
E(patient_graph)$weight <- E(patient_graph)$value

#Louvain method
louvain_communities <- cluster_louvain(graph = patient_graph, weights =E(patient_graph)$weight, resolution = 1.3)
louvain_communities$membership %>% table()
V(patient_graph)$louvain <-  louvain_communities$membership
com_table <- cbind.data.frame(id=V(patient_graph)$name, louvain=as.factor(louvain_communities$membership)) %>% as.data.frame()

rownames(com_table) <- com_table$id
com_table$id <- NULL
tail(com_table)

#visualize different order heatmap without hclust
com_table <- arrange(com_table, by=louvain)

louv_sim_matrix <- sim_matrix[rownames(com_table),rownames(com_table)]


#viz
colors <- rainbow(2)
color_pal <- viridis::viridis_pal(option = "B", direction = -1)(50)  
breaks <- seq(from = min(flat_matrix$value), to = max(flat_matrix$value), length.out = 11)
#legend <- c("474-1198","1198-5156","500-1000","1000-1800","> 1800")
legend_breaks <- breaks


palette.colors(palette = "Okabe-Ito")
colfunc <- colorRampPalette(c("palegreen", "seagreen4"))
colfunc(4)


ann_colors = list(
  gender = c(male="black", female="orange"),
  ajcc_pathologic_stage = c(stage1="#FFE4E1",stage2="#D898C4",stage3="#B14CA7",stage4="#8B008B"),
  ajcc_pathologic_t=c(T1="#ADD8E6",T2="#7390C7",T3="#3948A9",T4="#00008B",TX="white"),
  ajcc_pathologic_n=c(N0="#FFB6C1",N1="#D899A2",N2="#B17C83",N3="#8B5F65",NX="white"),
  ajcc_pathologic_m=c(M0="#98FB98",M1="#74D582",M1a="#51B06C",M1b="#2E8B57",MX="white"),
  louvain=c(`1`="#807dba",`2` = "#b2df8a",`3`="#33a02c", `4`="#e31a1c", `5`="#ff7f00", `6` = "#1f78b4"))



#png(filename = "10k_Patient_SSN_Sim_Heatmap.png", width = 480, height = 500, res = 120)
plot <- pheatmap(louv_sim_matrix, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames= F, na_col = "white", 
         fontsize = 7, color = color_pal,
         main = title, treeheight_col = 5, treeheight_row = 5 ,  annotation_col = clinic[,c("gender","ajcc_pathologic_stage", "ajcc_pathologic_t","ajcc_pathologic_n",
                                                                                            "ajcc_pathologic_m")], annotation_row = com_table,
         annotation_colors = ann_colors)
plot

# ggsave(
#   "10k_Patient_SSN_Sim_Heatmap_V2.png",
#   plot,
#   width = 8.5,
#   height = 8,
#   dpi = 150
# )
# 
# 
# 
# setwd("/datos/ot/plos")




#-------------------------------------------------------------------------------

#Extract cluster memberships for CoDiNA analysis and per-cluster consensus network building

cluster_tags <- levels(com_table[,1])
n_clusters <- length(cluster_tags)

#start loop per cluster:
for (i in 1:n_clusters) {
  
setwd(ssn_dir)
cluster_ids <- cluster_result[cluster_result[,1] == cluster_tags[i], ,drop=F] %>% rownames()
file_ids <- list.files(,pattern = "SSnet.tsv")
cluster_file_ids <- file_ids[sapply(file_ids, matches_pattern, patterns = cluster_ids)]

#format ids for CoDiNA
nr_samples <- length(cluster_file_ids) 
tcga_ids <- str_remove_all(string = cluster_file_ids,pattern = "_a_min_q-complete.tsv-SSnet.tsv")
tcga_ids <- str_remove_all(string = tcga_ids, pattern = "filtered-")
tcga_ids_mods <- gsub("-","_",tcga_ids)

#get nodes from expression matrix
node_list <- expr_matrix[,1]
empty_nodes <- cbind.data.frame(Var1=node_list,Var2=node_list,value=0)

#Create CoDiNA input list of networks
input <- c()

for (j in 1:nr_samples) {
  
  net <- fread(cluster_file_ids[j], sep = "\t") %>% as.data.frame()
  
  #comment next line if you want CoDiNA to run statistical analysis and filtering
  net$value <- ifelse(net$value >0, 1,-1)
  
  net <- rbind(net, empty_nodes)
  input <- c(input, list(net))
  
}

#run CoDiNA
diffnet <- MakeDiffNet(input, Code = tcga_ids_mods) %>% as.data.frame()

#Format diffnet object to output edge frequencies and arrange:

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

#Create codina outputs dir and save outputs
if (!file.exists(codina_outputs)) {
  # If it doesn't exist, create the directory
  dir.create(codina_outputs)
}

setwd(codina_outputs) 

diffnet <- diffnet[diffnet$row_counts >= quantile(diffnet$row_counts,.99),] #Get most conserved edges in a cluster
diffnet <- diffnet[,c("reg","tar","row_counts")]
output_cluster_name <- paste0("cluster-",as.character(cluster_tags[i]),"-diffnet.tsv")
#save
fwrite(diffnet,output_cluster_name, sep = "\t")
}


#Check Immune cell enrichment in clusters with Xcell
#-convert ensemble ids to Gene symbols

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- expr_matrix$genes

gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = gene_ids,
                      mart = ensembl)

gsymbol_expr_matr <- merge(expr_matrix, gene_symbols, by.x = "genes", by.y = "ensembl_gene_id", all.x = TRUE)
gsymbol_expr_matr <- gsymbol_expr_matr[!duplicated(gsymbol_expr_matr$external_gene_name),]
rownames(gsymbol_expr_matr) <- gsymbol_expr_matr$external_gene_name
gsymbol_expr_matr$genes <- NULL
gsymbol_expr_matr$external_gene_name <- NULL

#Run xCell 
xcell_obj <- xCellAnalysis(gsymbol_expr_matr)
xcell_obj <- xcell_obj[,rownames(com_table)]


#Perform kruskal wallis test on patient clusters and xCell outputs
t_xcell_obj <- t(xcell_obj) %>% as.data.frame() 

output <- cbind.data.frame(xcell_score = colnames(t_xcell_obj), p.val=NA)

for (i in 1:length(colnames(t_xcell_obj))) {
  
kruskal_table <- cooler_row_merge(com_table, t_xcell_obj[,i, drop=F])


kruskal_result <- kruskal.test(kruskal_table[,2] ~ kruskal_table[,1], data = kruskal_table)
output[i,2] <- kruskal_result$p.value

}

output$p.val <- p.adjust(output$p.val, method = "fdr")
output_significant <- output[output$p.val <= 0.001 ,] %>% arrange(., by=p.val)


#Posthoc Dunn test
kruskal_table_df <- cooler_row_merge(com_table, t_xcell_obj)

dunn_bxp <- \(fo, data, plot=TRUE, ...) {
  vars <- all.vars(fo)
  dnn <- rstatix::dunn_test(fo, data=data, p.adjust.method="BH")
  dnn_p <- rstatix::add_xy_position(dnn, x=vars[2])
  if (plot) {
    p <- ggpubr::ggboxplot(data, x=vars[2], y=vars[1], add = "jitter", color = "louvain",palette = c("#807dba","#b2df8a","#33a02c","#e31a1c","#ff7f00", "#1f78b4"), alpha = 0.5,label.rectangle = T) + 
      ggpubr::stat_pvalue_manual(dnn_p, label="p.adj.signif", hide.ns=TRUE, ...)
    print(p)
    #return(invisible(as.data.frame(dnn[1:8])))
    return(p)
  } else {
    return(as.data.frame(dnn[1:8]))
  }
}

#top immune
dunn_m1macro <- dunn_bxp(fo=`Macrophages M1` ~ louvain, data=kruskal_table_df)
dunn_adc <- dunn_bxp(fo=aDC ~ louvain, data=kruskal_table_df)
dunn_macro <- dunn_bxp(fo=Macrophages ~ louvain, data=kruskal_table_df)
dunn_th2 <- dunn_bxp(fo=`Th2 cells` ~ louvain, data=kruskal_table_df)
dunn_DC <- dunn_bxp(fo=DC ~ louvain, data=kruskal_table_df)
dunn_M2 <- dunn_bxp(fo=`Macrophages M2` ~ louvain, data=kruskal_table_df)
dunn_cDC <- dunn_bxp(fo=cDC ~ louvain, data=kruskal_table_df)
dunn_mono <- dunn_bxp(fo=Monocytes ~ louvain, data=kruskal_table_df)
dunn_iDC <- dunn_bxp(fo=iDC ~ louvain, data=kruskal_table_df)

#Immune, Stromal, HSC and Microenvironment score
dunn_Immune <- dunn_bxp(fo=ImmuneScore ~ louvain, data=kruskal_table_df)
dunn_Microenviron <- dunn_bxp(fo=MicroenvironmentScore ~ louvain, data=kruskal_table_df)
dunn_HSC <- dunn_bxp(fo=HSC ~ louvain, data=kruskal_table_df)
dunn_Stroma <-  dunn_bxp(fo=StromaScore ~ louvain, data=kruskal_table_df)

#top Stromal
dunn_adipocytes <- dunn_bxp(fo=Adipocytes ~ louvain, data=kruskal_table_df)
dunn_fibro <- dunn_bxp(fo=Fibroblasts ~ louvain, data=kruskal_table_df)
dunn_smooth <- dunn_bxp(fo=`Smooth muscle` ~ louvain, data = kruskal_table_df)
dunn_preadipo <- dunn_bxp(fo=Preadipocytes ~ louvain, data = kruskal_table_df)
dunn_myo <- dunn_bxp(fo=Myocytes ~ louvain, data = kruskal_table_df)
 

grid <- ggarrange(dunn_m1macro, dunn_adc, dunn_macro,dunn_th2,dunn_DC,dunn_adipocytes,dunn_fibro,dunn_smooth,dunn_preadipo,dunn_myo,common.legend = T,
          ncol = 5, nrow = 2,legend = T)

# grid
# ggsave(
#   "top5_immune_stromal.png",
#   grid,
#   width = 24,
#   height = 8,
#   dpi = 150
# )
# 


#Clinical features OVERREPRESENTATION TESTS------------------------------------------------------

enrich_cluster <- function(com_table,clinic, cat_variable){

colnames(com_table) <- "cluster"  
com_and_clinic <- cooler_row_merge(com_table, clinic)

n_factors <- com_and_clinic[[cat_variable]] %>% table() %>% length()
factors <- com_and_clinic[[cat_variable]] %>% as.factor() %>% levels() 

n_clusters <-  com_table$cluster %>% levels() %>% length()

p_vals = matrix(nrow = n_factors, ncol = n_clusters)


for (i in 1:n_clusters) {
  for (j in 1:n_factors) {
    
    #initialize mock variables
    aux_table <- com_and_clinic %>% mutate(., inCOM = "IN", inCAT= "IN")
    
    #Use Mock variables to fill in data necessary to build contingency table
    aux_table$inCOM = ifelse(aux_table$cluster == i, "IN","OUT")
    aux_table$inCAT = ifelse(aux_table[[cat_variable]] == factors[j],"IN","OUT")
    
    aux_table$inCOM <-  factor(aux_table$inCOM, levels=c("IN", "OUT") )
    aux_table$inCAT <-  factor(aux_table$inCAT, levels=c("IN", "OUT") )
    
    
    
    #Create contingency table and do fishers exact test
    aux_contin = aux_table %>% tabyl(.,inCAT,inCOM) 
    p_vals[j,i] = fisher.test(as.matrix(aux_contin[,-1]), alternative="greater")$p.value
    
    
  }
  
  
}

return(as.data.frame(p_vals))
}

stage_test <- enrich_cluster(com_table, clinic, "ajcc_pathologic_stage")
gender_test <- enrich_cluster(com_table, clinic, "gender")
vital_test <- enrich_cluster(com_table, clinic, "vital_status")
tissue_test <- enrich_cluster(com_table, clinic, "tissue_or_organ_of_origin")
tumor_test <- enrich_cluster(com_table, clinic, "ajcc_pathologic_t")
node_test <- enrich_cluster(com_table, clinic, "ajcc_pathologic_n")
met_test <- enrich_cluster(com_table, clinic, "ajcc_pathologic_m")


#Check survival differences between clusters
survival_table <- cooler_row_merge(com_table, clinic)
survival_table$vital_status <- ifelse(survival_table$vital_status=="Dead",1,0) 

lou.surv.plot_sim <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ louvain, data = survival_table) %>%
  ggsurvfit() +
  labs(
    x = "Days to last follow up",
    y = "Survival probability"
  ) + add_pvalue()


lou.surv.plot_sim

elify(gender_test) %>% p.adjust(.$value)

#Visualize Kruskal Wallis Results
Immune_cells <- c("aDC","B-cells","Basophils","CD4+ naive T−cells","CD4+ T−cells",
                  "CD4+ Tcm",
                  'CD4+ Tem',
                  'CD4+ memory T−cells',
                  'CD8+ naive T−cells',
                 'CD8+ T−cells',
                  'CD8+ Tcm',
                  'CD8+ Tem',
                  'cDC',
                  'Class−switched memory B−cells',
                  'DC',
                  'Eosinophils',
                  'iDC',
                  'Macrophages',
                  'Macrophages M1',
                  'Macrophages M2',
                  'Mast cells',
                  'Memory B−cells',
                  'Monocytes',
                  'naive B−cells',
                  'Neutrophils',
                  'NK cells',
                  'NKT',
                  'pDC',
                  'Plasma cells',
                  'pro B−cells',
                  'Tgd cells',
                  'Tregs',
                  'Th1 cells',
                  'Th2 cells')
Other_cells <- c('Astrocytes',
                'Epithelial cells',
                'Hepatocytes',
                'Keratinocytes',
               'Melanocytes',
                'Mesangial cells',
                'Neurons',
                'Sebocytes')
Stem_cells <- c('CLP',
                'CMP',
                'GMP',
                'HSC',
                'Megakaryocytes',
                'MPP',
                'Erythrocytes',
                'MEP',
                'Platelets')

Stromal_cells <- c('Adipocytes',
                   'Chondrocytes',
                   'Endothelial cells',
                   'Fibroblasts',
                   'ly Endothelial cells',
                   'MSC',
                   'mv Endothelial cells',
                   'Myocytes',
                   'Osteoblast',
                   'Pericytes',
                   'Preadipocytes',
                   'Skeletal muscle',
                   'Smooth muscle')
sig_immune <- output_significant %>% filter(., xcell_score %in% Immune_cells) %>% mutate(., Cell_type = "Immune") %>% .[order(.$p.val),]
sig_other <- output_significant %>% filter(., xcell_score %in% Other_cells) %>% mutate(., Cell_type = "Other")  %>% .[order(-.$p.val),]
sig_stem <- output_significant %>% filter(., xcell_score %in% Stem_cells) %>% mutate(., Cell_type = "Stem")  %>% .[order(-.$p.val),]
sig_stromal <- output_significant %>% filter(., xcell_score %in% Stromal_cells) %>% mutate(., Cell_type = "Stromal")  %>% .[order(.$p.val),]

significant_w_type <- rbind.data.frame(sig_immune,sig_stromal) %>% mutate(.,  `-log(FDR)`= -log(.$p.val))

dot_chart <- ggdotchart(significant_w_type, x = "xcell_score", y = "`-log(FDR)`",
           color = "Cell_type",                                
           palette = c("#4E79A7" ,"#F28E2B" ),                      
           rotate = TRUE,                                
           dot.size = 5,                                 
           y.text.col = TRUE,
           sorting = "none",
           ggtheme = theme_pubr()                        
          ) +
  theme_cleveland()                                      



dot_chart

kuskal_posthoc <- ggarrange(dot_chart, dunn_Immune, dunn_Stroma,dunn_Microenviron, common.legend = F,
          ncol = 2, nrow = 2,legend = T, labels = c("A","B","C","D"))


# ggsave(
#   "Cell_types_Kruskal_grid_v2.png",
#   kuskal_posthoc,
#   width = 10,
#   height = 10,
#   dpi = 150
# )



# ggsave(
#   "Cell_types_Kruskal.png",
#   dot_chart,
#   width = 7,
#   height = 8.5,
#   dpi = 150
# )

#Plot heatmap of median xcell scores

kruskal_table_df %>% filter(., louvain == "1")
  
cls1_cell_medians <-  kruskal_table_df %>% filter(., louvain=="1") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()
cls2_cell_medians <-  kruskal_table_df %>% filter(., louvain=="2") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()
cls3_cell_medians <-  kruskal_table_df %>% filter(., louvain=="3") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()
cls4_cell_medians <-  kruskal_table_df %>% filter(., louvain=="4") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()
cls5_cell_medians <-  kruskal_table_df %>% filter(., louvain=="5") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()
cls6_cell_medians <-  kruskal_table_df %>% filter(., louvain=="6") %>% .[,-1] %>% apply(., MARGIN = 2,FUN = median) %>% as.data.frame()

merged_cls_medians <- data.frame(Clsuter_1 = cls1_cell_medians,
                                       Clsuter_2 = cls2_cell_medians,
                                       Clsuter_3 = cls3_cell_medians,
                                       Clsuter_4 = cls4_cell_medians,
                                       Clsuter_5 = cls5_cell_medians,
                                       Clsuter_6 = cls6_cell_medians)

colnames(merged_cls_medians) <- c("Cluster_1","Cluster_2", "Cluster_3","Cluster_4", "Cluster_5", "Cluster_6")
merged_cls_medians <- merged_cls_medians %>% .[rownames(.) %in% significant_w_type$xcell_score,] %>% .[significant_w_type$xcell_score,] %>% apply(., MARGIN =1, FUN=CoDiNA::normalize) %>% t()
  
all_cells_sigf <- pheatmap(merged_cls_medians, scale = "row", cluster_cols = F, cluster_rows = T)

# ggsave(
#   "heatmap_significant_celltypes.png",
#   all_cells_sigf,
#   width = 5,
#   height = 7,
#   dpi = 150
# )
# 

#overlap plot with sankey diagram

ssn_dir_10k <- "/datos/ot/plos/aracne-multicore/launch/filtered_10k_SSnets_surv_luadtp"
ssn_dir_50k <- "/datos/ot/plos/aracne-multicore/launch/filtered_50k_SSnets_surv_luadtp"
ssn_dir_100k <- "/datos/ot/plos/aracne-multicore/launch/filtered_100k_SSnets_surv_luadtp"


setwd(dir)

#load similarity matrices and format
sim_matrix_10 <- fread("filtered10k_surv_SSnets_similarity.tsv", sep = "\t") %>% as.data.frame()
sim_matrix_50 <- fread("filtered50k_surv_SSnets_similarity.tsv", sep = "\t") %>% as.data.frame()
sim_matrix_100 <- fread("filtered100k_surv_SSnets_similarity.tsv", sep = "\t") %>% as.data.frame()

sample_id <- colnames(expr_matrix)
faulty <- sample_id[2]

rownames(sim_matrix_10) <- colnames(sim_matrix_10)
rownames(sim_matrix_50) <- colnames(sim_matrix_50)
rownames(sim_matrix_100) <- colnames(sim_matrix_100)

sim_matrix_10 <- sim_matrix_10[! rownames(sim_matrix_10) == faulty, ! colnames(sim_matrix_10) == faulty]
sim_matrix_50 <- sim_matrix_50[! rownames(sim_matrix_50) == faulty, ! colnames(sim_matrix_50) == faulty]
sim_matrix_100 <- sim_matrix_100[! rownames(sim_matrix_100) == faulty, ! colnames(sim_matrix_100) == faulty]

file_ids <- colnames(sim_matrix_10)
tcga_ids <- str_remove_all(string = file_ids, 
                           pattern = "_a_min_q-complete.tsv-SSnet.tsv")
tcga_ids <- str_remove_all(string = tcga_ids, 
                           pattern = "filtered-")

colnames(sim_matrix_10) <- tcga_ids
rownames(sim_matrix_10) <- tcga_ids

colnames(sim_matrix_50) <- tcga_ids
rownames(sim_matrix_50) <- tcga_ids

colnames(sim_matrix_100) <- tcga_ids
rownames(sim_matrix_100) <- tcga_ids

#transform similarity matrix form # of edges to jaccard similarity index
n_edges_10 <- 10000
n_edges_50 <- 50000
n_edges_100 <- 100000

jaccard_mat_sim <- function(sim_matrix, n_edges){
x <- sim_matrix/((n_edges*2)-sim_matrix)
return(x)
}

jaccard_10 <- jaccard_mat_sim(sim_matrix_10,  n_edges_10)
jaccard_50 <- jaccard_mat_sim(sim_matrix_50,  n_edges_50)
jaccard_100 <- jaccard_mat_sim(sim_matrix_100,  n_edges_100)

flat_matrix_10 <- elify(as.matrix(jaccard_10)) %>% .[order(- .[,3]),] 
flat_matrix_50 <- elify(as.matrix(jaccard_50)) %>% .[order(- .[,3]),] 
flat_matrix_100 <- elify(as.matrix(jaccard_100)) %>% .[order(- .[,3]),] 

graph_clust <- function(flat_matrix){
patient_graph <- graph_from_data_frame(flat_matrix, directed = F)
E(patient_graph)$weight <- E(patient_graph)$value

#Louvain method
louvain_communities <- cluster_louvain(graph = patient_graph, weights =E(patient_graph)$weight, resolution = 1.3)
louvain_communities$membership %>% table()
V(patient_graph)$louvain <-  louvain_communities$membership
com_table <- cbind.data.frame(id=V(patient_graph)$name, louvain=as.factor(louvain_communities$membership)) %>% as.data.frame()

rownames(com_table) <- com_table$id
com_table$id <- NULL
return(com_table)
}

louv_clust_10 <- graph_clust(flat_matrix_10)
louv_clust_50 <- graph_clust(flat_matrix_50)
louv_clust_100 <- graph_clust(flat_matrix_100)

merged_coms <- cooler_row_merge(louv_clust_10, louv_clust_50) %>% cooler_row_merge(., louv_clust_100) 
colnames(merged_coms) <- c("10 000 edges clustering", "50 000 edges clustering", "100 000 edges clustering")


df <- merged_coms %>%
  make_long(`10 000 edges clustering`, `50 000 edges clustering`, `100 000 edges clustering`)
df



dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()


df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

pl <- ggplot(df2, aes(x = x
                      , next_x = next_x
                      , node = node
                      , next_node = next_node
                      , fill = factor(node)
                      
                      , label = paste0("Cluster ",node)
)
) 
pl <- pl +geom_sankey(flow.alpha = 0.35,color="grey40",show.legend = TRUE)
pl <- pl +geom_sankey_label(size = 5, color = "white", fill= "gray40", hjust = -0.2)
pl <- pl + theme(legend.position = "none")
pl <- pl + theme_bw()
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(fill = 'Nodes')


pl

ggsave(
  "sankey_overlap.png",
  pl,
  width = 15,
  height = 10,
  dpi = 250
)



