#!/bin/bash

DIR=$1

COUNT_REPORT="${DIR}/profileTracker.txt"

echo -e "ReadsIn\tDiscordantOverlappingReadFragment\tDisparateCoord\tHasMS\tLowCoverage\tLowMapQ\tMsInPCRDuplicates\tNonReference\tOverlappingReadFragment\tReported\tUnmapped\tWrongUnit\tReadsOut" > $COUNT_REPORT
for fam in $(find $DIR -maxdepth 1 -name "auSSC*" -type d)
do
	for person in $(find $fam -name "SSC*" -type d)
	do
		PERSON_ID=$(basename $person)
		echo -ne "$PERSON_ID" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep "reads in" {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep discordantOverlappingReadFragment {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep disparateCoord {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep hasMS {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep lowCoverage {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep lowMapQ {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep msInPCRDuplicates {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep nonReference {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep overlappingReadFragment {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep reported {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep unmapped {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -ne "\t$(find $person -name *.profile.count -exec grep wrongUnit {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
		echo -e "\t$(find $person -name *.profile.count -exec grep out {} \; | cut -f2 -d ':' | sed 's/ //g')" >> $COUNT_REPORT
	done
done

echo "Count tracker can be found at $COUNT_REPORT"
