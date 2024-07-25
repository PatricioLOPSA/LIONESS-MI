library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

#INPUT
projname <- "TCGA-LUAD"

#PROCESS
query <- GDCquery(project = projname,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query)
rnas <- GDCprepare(query = query, summarizedExperiment = TRUE)

file_name <- paste(projname,"rnas-tcga.rds", sep = "-" )
saveRDS(rnas, file = file_name)