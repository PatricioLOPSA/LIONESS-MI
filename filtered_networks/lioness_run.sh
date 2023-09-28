# Batch run
# for i in $(ls *.tsv); do bash run.sh $i; done > salida &

partools="$(pwd)/../parallel"
aracnebin="../bin/aracne2"

[[ ! -f $aracnebin ]] \
	&& echo "No binary: $aracnebin" && exit 15

[[ ! -d $partools ]] \
	&& echo "No parallel tools: ../bin/aracne2" && exit 15

ftsv=$1

[[ $ftsv == "" ]] && echo "need .tsv file" && exit 15

# ftsv="norm-Control.tsv" 
echo "Processing MI calculations for: $ftsv ..."
nom=$(echo $ftsv | cut -d. -f 1)

awk '{print $1}' $ftsv > node.list
cname=$(head -1 node.list)
echo "Column Index Name: $cname"

SECONDS=0
python ${partools}/aracne-par.py $ftsv node.list $cname $(nproc) &> aracne.log 
echo "ARACNe time: $(echo $SECONDS/60 | bc -l) minutes."

SECONDS=0
n=$( (cd adj; ls) | head -1 | cut -d'.' -f 2 )
echo "Parameters to join: $nom $n node.list $cname"
python ${partools}/joinadj.py $nom $n node.list $cname
echo "join ADJ matriz time: $(echo $SECONDS/60 | bc -l) minutes."

echo "Moving adjancy matrix"
mv adj/mat.adj .

SECONDS=0
python ${partools}/adj2tsv.py mat.adj
mv mat-complete.tsv  ${nom}-complete.tsv
echo "Creating Complete MAtrix: $(echo $SECONDS/60 | bc -l) minutes."


#SECONDS=0
#python ${partools}/adj2sif.py > ${nom}.sif
#echo "Creating SIF: $(echo $SECONDS/60 | bc -l) minutes."

# SECONDS=0
# # sort -r -k3,3 ${nom}.sif | head -10000 > ${nom}-IM-1e5.txt
# sort -r -k3,3 ${nom}.sif > ${nom}.sort
# echo "Sorting: $(echo $SECONDS/60 | bc -l) minutes."
# # for i in $(\ls *.sort); do n=${i%.*}; n=${n/matexp/IM}; head -10000 $i > "$n-10k.txt"; done

rm -rf adj log mat.adj node.list 
