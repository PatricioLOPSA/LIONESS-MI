# edited from https://github.com/josemaz/kirc-mirna/blob/main/R/biomart.R

require(biomaRt)
require(dplyr)
require(tidyverse)

ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
print("Downloading BIOMART data")

bm <- getBM(attributes=c('ensembl_gene_id',
                         'description',
                         'gene_biotype',
                         'percentage_gene_gc_content',
                         'start_position',
                         'end_position',
                         'transcript_length',
                         'chromosome_name',
                         'band',
                         'external_gene_name'),
            mart = ensembl,
            verbose = FALSE)
chroms <- c(1:22,"X","Y")
dtmp <- subset(bm, chromosome_name %in% chroms)
dtmp <- subset(dtmp, gene_biotype %in% 
                 c("protein_coding"))
dtmp <- dtmp[!duplicated(dtmp$ensembl_gene_id), ]
dtmp$geneLength <- dtmp$end_position - dtmp$start_position
dtmp <- rename(dtmp, gcContent = percentage_gene_gc_content)
dtmp <- rename(dtmp, chr = chromosome_name)
dtmp <- rename(dtmp, gene_name = external_gene_name)
dtmp <- dtmp  %>% drop_na(gene_name)
dtmp <- dtmp[dtmp$gene_name != "", ]

write.csv(dtmp, file = "biomart.csv", row.names = FALSE)
