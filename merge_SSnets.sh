#!/bin/bash


dir=$1

cd $dir

for file in $(ls *SSnet.tsv); do


	id=$(echo $file | cut -d_ -f 1)
	valuefile=$(echo "$id"_aux.tsv)

	cut -f3 $file > $valuefile
    
	sed 1s/value/$id/ $valuefile > ${id}_holder.tsv
	rm $valuefile
	
done

paste <(cut -f 1,2 $(ls *SSnet.tsv | head -1)) *holder* > merged_lioness_nets.tsv
rm *holder*


