#!/bin/bash

DIR=$1

COUNT_REPORT="${DIR}/replicateTracker.txt"

[ -f $COUNT_REPORT ] && rm $COUNT_REPORT

echo -e "person\treplicate.sets\tdiscordant.replicate.sets\tnum.pairs" > $COUNT_REPORT
for fam in $(find $DIR -maxdepth 1 -name "auSSC*" -type d)
do
	for person in $(find $fam -name "SSC*" -type d)
	do
		PERSON_ID=$(basename $person)
		DUP_FILE=$(find $person -name "*.dups.txt.gz")
		TOTAL=$(gunzip -c $DUP_FILE | wc -l)
		DISC=$(gunzip -c $DUP_FILE | awk '$8 > 0 {c++};END{print c}')
		NUMPAIR=$(gunzip -c $DUP_FILE | awk '$9 == 2 {c++};END{print c}')
		echo -e "$PERSON_ID\t$TOTAL\t$DISC\t$NUMPAIR" >> $COUNT_REPORT
	done
done

echo "PCR replicate tracker can be found at $COUNT_REPORT"
