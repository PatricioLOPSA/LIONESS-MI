# LIONESS-MI Networks
This project contains the tools needed to reconstruct single sample networks using  a parallelized implementation of ARACNE (https://github.com/josemaz/aracne-multicore) and the LIONESS equation. 

- ``lioness_prepare.R`` Takes a complete ($\alpha$) expression matrix  as input, and iteratively computes and stores all $\alpha - q$ matrices.  First column of input matrix must be a list of gene ids and should be stored as a .tsv file


- ``lioness_run.sh`` Calls the multicore implementation of ARACNe2. Takes an expression matrix and computes de Mutual Information (MI) for all gene pairs. The output is an adjacency matrix stored as .tsv file.

- ``lioness_networks.R`` Computes the LIONESS equation for any $\alpha - q$ matrix. The resulting network is filtered to keep the top 100k edges filtered by the absolute value of their LIONESS' score.

- ``lioness.sh`` Calls all the above mentioned scripts to iteratively calculate all single sample networks using ARACNe2 .

- ``CoDiNA.R`` Uses the CoDiNA R package to join all single sample networks, and search for patient specific edges. 
