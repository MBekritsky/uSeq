#!/bin/bash

DIR=$1
NUM_LOCI=$2
DN_SCRIPT_PATH="/mnt/wigclust10/data/safe/bekritsk/simons/scripts/genotyper/emGenotyper/loadAndGraphDenovos.R"

[ -z $NUM_LOCI ] && NUM_LOCI=10000

TOTAL_LOCI=$(awk '$8 == 0 && $1 ~ /[0-9]$/' "${DIR}/allele_matrix_info.txt" | wc -l)
SGEERR="${DIR}/DnGraphSGEERR/"
SGEOUT="${DIR}/DnGraphSGEOUT/"

[ ! -d $SGEERR ] && $(mkdir $SGEERR)
[ ! -d $SGEOUT ] && $(mkdir $SGEOUT)

for i in $(seq 1 $NUM_LOCI $TOTAL_LOCI)
do
ERRFILE="${SGEERR}dn${i}.err"
OUTFILE="${SGEOUT}dn${i}.out"
OUTPUT=$(qsub -q all.q@wigclust10.cshl.edu -S /bin/bash -e $ERRFILE -o $OUTFILE $DN_SCRIPT_PATH --directory=$DIR --start.locus=$i --num.loci=$NUM_LOCI --script.file=$DN_SCRIPT_PATH)
echo "Submitted loci $i to $(($i + $NUM_LOCI))"
echo "$OUTPUT"
done
