# remove previous outputs
rm *_a_min_q-complete.tsv *_a_min_q.tsv


exp=$1
output=$2

nom=$(echo $exp | cut -d. -f 1)

Rscript lioness_prepare.R $exp

bash lioness_run.sh $exp

mv ${nom}-complete.tsv alpha_net.tsv

for i in $(ls *_a_min_q.tsv); do bash lioness_run.sh $i; done

for i in $(ls *_a_min_q-complete.tsv); do Rscript lioness_networks_posints.R $exp $i ; done

rm *_a_min_q-complete.tsv *_a_min_q.tsv

mkdir $output
mv *alpha_net.tsv *100k-SSnet.tsv $output
