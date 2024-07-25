require(SummarizedExperiment)
require(dplyr)
require(stringr)


luad <- readRDS("TCGA-LUAD-rnas-tcga.rds")

#setwd(/datos/ot/plos)
bm <- read.csv("biomart.csv")





###################################################################3
# source: https://github.com/josemaz/kirc-mirna/blob/main/R/tools.R

# FILTER
#! Filter by minimum mean and number of zeros
filtro.zeros.means <- function(m1) {
  print("Input (dim): ")
  print(dim(m1))
  threshold <- round(dim(m1)[2]/2) #Zero-prevalence filter
  print(paste0("Threshold: ",threshold))
  m1 <- m1[rowSums(m1 == 0) <= threshold, ] 
  #no pasen el umbral
  
  print(paste0("Rows After Zeros (dim): ",dim(m1)[1]))
  m1 <- m1[rowMeans(m1) >= 10, ]    # minimum mean counts filtering
  print(paste0("Rows After means (dim): ",dim(m1)[1]))
  return(m1)
}

####################################################################

####################################################################
#!-- RNA Samples Cleaning
clean.rna <- function(RNA.exp,dat.bm){
  # TODO: print correctly dims reports
  #!-- Filtros
  print(dim(RNA.exp))
  lstg <- c("Primary solid Tumor","Solid Tissue Normal")
  rnaseq <- RNA.exp[,RNA.exp$definition %in% lstg]
  dim(rnaseq)
  lstg <- c("Stage IA","Stage I","Stage II","Stage III",
            "Stage IB","Stage IIA","Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IV")
  # TODO: Rename "stage i" -> stage1, etc.
  # rnaseq <- rnaseq[,rnaseq$tumor_stage %in% lstg]
  #lstg <- c("Stage I","Stage II","Stage III","Stage IV")
  rnaseq <- rnaseq[,rnaseq$ajcc_pathologic_stage %in% lstg]
  print(dim(rnaseq))
  
  #!-- Mutate and Changes
  clinic <- colData(rnaseq)
  rnaseq$grupo <- rnaseq$ajcc_pathologic_stage
  rnaseq$grupo[clinic$sample_type_id=="11"] <- "ctrl"
  rnaseq$grupo[rnaseq$grupo=="Stage I"] <- "stage1"
  rnaseq$grupo[rnaseq$grupo=="Stage IA"] <- "stage1"
  rnaseq$grupo[rnaseq$grupo=="Stage IB"] <- "stage1"
  rnaseq$grupo[rnaseq$grupo=="Stage II"] <- "stage2"
  rnaseq$grupo[rnaseq$grupo=="Stage IIA"] <- "stage2"
  rnaseq$grupo[rnaseq$grupo=="Stage IIB"] <- "stage2"
  rnaseq$grupo[rnaseq$grupo=="Stage III"] <- "stage3_4"
  rnaseq$grupo[rnaseq$grupo=="Stage IIIA"] <- "stage3_4"
  rnaseq$grupo[rnaseq$grupo=="Stage IIIB"] <- "stage3_4"
  rnaseq$grupo[rnaseq$grupo=="Stage IV"] <- "stage3_4"
  rnaseq$grupo <- as.factor(rnaseq$grupo)
  
  
  #--! Eliminate duplicates by group
  print(dim(rnaseq))
  print(paste0("Sample dups: ",ncol(rnaseq[,duplicated(rnaseq$sample)])))
  clinic <- colData(rnaseq)
  levels(clinic$grupo)
  d1 <- clinic[clinic$definition == "Primary solid Tumor",]
  d1 <- d1[!duplicated(d1$patient),] 
  d2 <- clinic[clinic$definition == "Solid Tissue Normal",]
  d2 <- d2[!duplicated(d2$patient),] 
  rnaseq <- rnaseq[,c(rownames(d1),rownames(d2))]
  rm(d1,d2,lstg,clinic)
  print(dim(rnaseq))
  
  #!-- GENES Cleaning
  m <- filtro.zeros.means(assay(rnaseq))
  rnaseq <- rnaseq[rownames(m),]
  dim(rnaseq)
  
  
  #!-- Annotation
  # bm <- read.csv("pipeline/biomart.csv")
  bm <- dat.bm[dat.bm$gene_biotype == "protein_coding",]
  bm <- bm[!duplicated(bm$ensembl_gene_id),]
  ensmbl_ids_version <- rownames(rnaseq) #-------PLS
  rownames(rnaseq) <- substring(ensmbl_ids_version,1,15)#-----PLS
  rnaseq <- rnaseq[rownames(rnaseq) %in% bm$ensembl_gene_id, ]#------PLS
  bm <- bm[bm$ensembl_gene_id %in% rownames(rnaseq),]
  rowData(rnaseq)$ensembl_gene_id <- rownames(rnaseq)
  rowData(rnaseq) <- merge(rowData(rnaseq),bm,by="ensembl_gene_id")
  print(dim(rnaseq))
  return(rnaseq)
}


###################################################################3

luad_clean_stg <- clean.rna(luad,bm) #(16289,575)



#setwd("C:/Users/patri/data/RDS")

saveRDS(luad_clean_stg, "rnaseq-luad-clean.rds")

