#!/bin/bash

DIR=$1

for fam in $(find $DIR -maxdepth 1 -name "auSSC*" -type d)
do
	FAMID=$(basename $fam)
	for person in $(find $fam -maxdepth 1 -name "SSC*" -type d)
	do
		PERSONID=$(basename $person)
		LOGPATH="${person}/logs/*"
		SGEERRPATH="${person}/SGEERR/*"
		SUBMITERR="${person}/${PERSONID}.err"
		ERRORS=$(grep "rror" $LOGPATH; grep "rror" $SGEERRPATH; grep "rror" $SUBMITERR)

		echo "$FAMID $PERSONID"
		if [ -n "$ERRORS" ]	
		then
			echo "ERRORS"
			echo "$ERRORS"
		fi
		
		if [ -z "$ERRORS" ]
		then
			echo "No errors"
		fi

	done
done
