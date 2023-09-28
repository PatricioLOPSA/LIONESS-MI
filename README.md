# LIONESS-MI Networks
This project contains the tools needed to reconstruct single sample networks using a parallelized implementation of ARACNE (https://github.com/josemaz/aracne-multicore) and the LIONESS equation. 

- ``get_aminq_mat.R`` Takes a complete ($\alpha$) expression matrix  as input, and iteratively removes sample $q$ and stores all $\alpha - q$ matrices.  First column of input matrix must be a list of gene ids and should be stored as a .tsv file



- ``aracne_run.sh`` Calls the multicore implementation of ARACNe2. Takes an expression matrix and desired number of processors as arguments and computes the Mutual Information (MI) for all gene pairs. The output is an adjacency matrix stored as .tsv file.



- ``get_SS_net.R``  Uses the LIONESS equation for any $\alpha - q$ matrix. The resulting networks are a list of all possible edges and their LIONESS score filtered to keep the upper triangular part of a symmetric adjacency matrix. All networks are stored as separate files for downstream analyses.



- ``lioness_aracne.sh`` Calls all the above mentioned scripts to calculate all single sample networks using ARACNe2.



- ``CoDiNA.R`` Uses the CoDiNA R package to join all single sample networks, and search for patient specific edges (see scripts in ``filtered_networks``. 

## Hands on

1. Follow the instructions in https://github.com/josemaz/aracne-multicore to install ARACNe-multicore.

2.  Clone this repository, and move all files and expression matrix of interest into the  ``launch`` directory.

3.  Run ``lioness_aracne.sh`` with the corresponding arguments:
- first argument must be the expression matrix of interest that meets the above mentioned criteria.
- second argument is the desired name of the output directory that will store all sample specific and aggregate networks.
- third argument is the desired number of processors that ARACNE should use in the network reconstruction process.
```bash
bash lioness_aracne.sh expression_matrix.tsv outputdir 20
```
#### Optional

4. Run Co-expression Differential Network Analysis (CoDiNA) R script to merge all networks into a single data frame and look for patient specific edges. 
```bash
Rscript CoDiNA.R outputdir expression_matrix.tsv
```
CoDiNA output is stored as ``outputdir-CoDiNA-diffnet.tsv``.


