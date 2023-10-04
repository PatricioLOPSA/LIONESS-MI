# remove previous outputs
rm *_a_min_q.tsv


exp=$1
output=$2
nprocs=$3


nom=$(echo $exp | cut -d. -f 1)

Rscript get_aminq_mat.R $exp

bash aracne_run.sh $exp $nprocs

mv ${nom}-complete.tsv alpha_net.tsv

for i in $(ls *_a_min_q.tsv); do if [ ! -e $(echo $i | cut -d. -f 1)-complete.tsv ]; then bash aracne_run.sh $i $nprocs else echo "a_min_q net already exists"; fi; done

for i in $(ls *_a_min_q-complete.tsv); do Rscript get_SS_net.R $exp $i ; done

rm *_a_min_q-complete.tsv *_a_min_q.tsv

mkdir $output
mv *alpha_net.tsv *-SSnet.tsv $output
