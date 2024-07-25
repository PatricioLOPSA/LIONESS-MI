
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

setwd("/datos/ot/plos/")
data <- readRDS("rnaseq-luad-clinicstg.rds")

#Work with primary solid tumor only
tp <- data[,data$grupo == "Primary solid Tumor"]

#remove patients with NA as follow up data
tp <- tp[,!is.na(tp$days_to_last_follow_up)]

#Left censor patients who were not followed up for more than a year or did not experience event
surv_threshold = 365


censored <- tp[,(tp$days_to_last_follow_up < surv_threshold) & (tp$vital_status == "Alive")]
tp_filt_surv <- tp[,!(colnames(tp) %in% colnames(censored))]

surv_luadtp_matrix <- assay(tp_filt_surv) %>% as.data.frame()
surv_luadtp_matrix <- cbind.data.frame(genes = rownames(surv_luadtp_matrix), surv_luadtp_matrix)
surv_luadtp_matrix[1:5,1:5]

fwrite(surv_luadtp_matrix, "survival_luadtp_matrix.tsv", sep = "\t", row.names = F)

