#!/bin/bash

DIR=$1

COUNT_REPORT="${DIR}/inOutTracker.txt"

echo -e "person\tscan.in\tscan.out\treindex.in\treindex.out\tprofile.in" > $COUNT_REPORT
for fam in $(find $DIR -maxdepth 1 -name "auSSC*" -type d)
do
	for person in $(find $fam -name "SSC*" -type d)
	do
		PERSON_ID=$(basename $person)
		SCAN_IN=$(find $person -name *.scan.count -exec grep in {} \; | awk '{c+=$3};END{print c}')
		SCAN_OUT=$(find $person -name *.scan.count -exec grep out {} \; | awk '{c+=$3};END{print c}')
		REINDEX_IN=$(find $person -name *reindex.count -exec grep in {} \; | awk '{c+=$3};END{print c}')
		REINDEX_OUT=$(find $person -name *reindex.count -exec grep out {} \; | awk '{c+=$3};END{print c}')
		PROFILE_IN=$(find $person -name *profile.count -exec grep in {} \; | awk '{c+=$3};END{print c}')
		echo -e "$PERSON_ID\t$SCAN_IN\t$SCAN_OUT\t$REINDEX_IN\t$REINDEX_OUT\t$PROFILE_IN" >> $COUNT_REPORT
	done
done

echo "Count tracker can be found at $COUNT_REPORT"
