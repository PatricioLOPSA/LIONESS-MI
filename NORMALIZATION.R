require(SummarizedExperiment)
require(NOISeq)
require(EDASeq)
require(DESeq2)
require(dplyr)


luad <- readRDS("rnaseq-luad-clean.rds")


###############################################################################
# source: https://github.com/josemaz/kirc-mirna/blob/main/R/tools.R
#!-- PCA plotting by group
pca.grupo <- function(dat, fout = NULL){
  d2 <- data.frame(tejido = colData(dat)$grupo)
  rownames(d2) <- colnames(dat)
  mydat = NOISeq::readData( assay(dat) , factors = d2)
  myPCA = dat(mydat, type = "PCA", logtransf = F)
  if(!is.null(fout)){
    print(paste0("Writing in: ",fout))
    png(fout)
  }
  explo.plot(myPCA, factor = "tejido", plottype = "scores")
  if(!is.null(fout))dev.off()
}

###############################################################################

###############################################################################
# Normalization and Bias correct of RNAseq
norm.rna <- function(dat.rna, dout){
  pca.grupo(dat.rna,fout = paste0(dout,"/PCA-rna-BeforeNorm.png"))
  fac <- data.frame(tejido=dat.rna$grupo, 
                    row.names=colnames(dat.rna))
  #! Pre Normalization
  ln.data <- withinLaneNormalization(assay(dat.rna), 
                                     rowData(dat.rna)$geneLength, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , rowData(dat.rna)$gcContent,
                                      which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData( norm.counts, factors = fac)
  #! Post Normalization
  assay(dat.rna) <- exprs(noiseqData)
  pca.grupo(dat.rna,fout = paste0(dout,"/PCA-rna-AfterNorm.png"))
  return(dat.rna)
}

###############################################################################


dout <- "plots"
luadnorm <- norm.rna(luad,dout)


saveRDS(luadnorm, "rnaseq-luad-clinicstg.rds")


# validated by PhD. JE
