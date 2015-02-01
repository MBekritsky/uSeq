#!/bin/bash

DIR=$1
NUM_LOCI=$2
CG_SCRIPT_PATH="/mnt/wigclust4/home/bekritsk/tools/microsatellite/emGenotyper/chunkGenotyper.R"

[ -z $NUM_LOCI ] && NUM_LOCI=10000
TOTAL_LOCI=$(awk '$8 == 0 && $1 ~ /[0-9]$/' "${DIR}/allele_matrix_info.txt" | wc -l)

SGEERR="${DIR}/ChunkGenotypeSGEERR/"
SGEOUT="${DIR}/ChunkGenotypeSGEOUT/"

[ ! -d $SGEERR ] && $(mkdir $SGEERR)
[ ! -d $SGEOUT ] && $(mkdir $SGEOUT)

for i in $(seq 1 $NUM_LOCI $TOTAL_LOCI)
do
ERRFILE="${SGEERR}cg${i}.err"
OUTFILE="${SGEOUT}cg${i}.out"
OUTPUT=$(qsub -S /bin/bash -e $ERRFILE -o $OUTFILE $CG_SCRIPT_PATH --directory=$DIR --start.locus=$i --num.loci=$NUM_LOCI --script.file=$CG_SCRIPT_PATH)
echo "Submitted loci $i to $(($i + $NUM_LOCI))"
echo "$OUTPUT"
done
