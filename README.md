# LIONESS-MI Networks
This project contains the tools needed to reconstruct single sample networks using a parallelized implementation of ARACNE (https://github.com/josemaz/aracne-multicore) and the LIONESS equation. 

- ``lioness_prepare.R`` Takes a complete ($\alpha$) expression matrix  as input, and iteratively removes sample $q$ and stores all $\alpha - q$ matrices.  First column of input matrix must be a list of gene ids and should be stored as a .tsv file



- ``lioness_run.sh`` Calls the multicore implementation of ARACNe2. Takes an expression matrix and computes the Mutual Information (MI) for all gene pairs. The output is an adjacency matrix stored as .tsv file.



- ``lioness_networks.R``  Uses the LIONESS equation for any $\alpha - q$ matrix. The resulting network is filtered to keep the top 100k edges based on the absolute value of their LIONESS' score.



- ``lioness.sh`` Calls all the above mentioned scripts to calculate all single sample networks using ARACNe2.



- ``CoDiNA.R`` Uses the CoDiNA R package to join all single sample networks, and search for patient specific edges. 

## Hands on

1. Follow the instructions in https://github.com/josemaz/aracne-multicore to install ARACNe-multicore.

2.  Clone this repository, and move all files and expression matrix of interest into the  ``launch`` directory.

3.  Run ``lioness.sh`` with the corresponding arguments:
- first argument must be the expression matrix of interest that meets the above mentioned criteria.
- second argument is the desired name of the output directory that will store all sample specific and aggregate networks.
```bash
bash lioness.sh expression_matrix.tsv outputdir
```
#### Optional

4. Run Co-expression Differential Network Analysis (CoDiNA) R script to merge all networks into a single data frame and look for patient specific edges. 
```bash
Rscript CoDiNA.R outputdir expression_matrix.tsv
```
CoDiNA output is stored as ``outputdir-CoDiNA-diffnet.tsv``.


