#!/bin/bash

INFILE=$1
DIR=$2

[ -z "$DIR" ] && DIR="."

REPORT_FILE="${DIR}/copyReport_$(date +%d%m%Y).txt"
[ -f $REPORT_FILE ] && rm $REPORT_FILE 


sed 1d $INFILE | while read ENTRY #skips the header line in the report file
do
	if [[ $(echo "$ENTRY" | cut -f30 | grep -c "joe") > 0  ]]
	then
		echo "Skipping joe" >> $REPORT_FILE
		continue
	fi
	
	if [[ $(echo "$ENTRY" | cut -f30 | grep -c "someone") > 0  ]]
	then
		echo "Skipping \"someone\"" >> $REPORT_FILE
		continue
	fi

	if [[ $(echo "$ENTRY" | cut -f26 | grep -c "NULL") > 0  ]]
	then
		echo "Skipping NULL family" >> $REPORT_FILE
		continue
	fi
		
	SEQINFO=$(echo "$ENTRY" | cut -f2)
	FAMILY_ID=$(echo "$ENTRY" | cut -f26)
	REL=$(echo "$ENTRY" | cut -f28)
	INDIVIDUAL_ID=$(echo "$ENTRY" | cut -f30)

	FC=$(echo $SEQINFO | cut -f1 -d '-')
	LANE=$(echo $SEQINFO | cut -f2 -d '-')
	BC=$(echo "$ENTRY" | cut -f3)

	FILEPATH=$(echo "$ENTRY" | cut -f13)
	FILEPATH="${FILEPATH}/bc${BC}"

	INFO="${FC}${LANE}"

	LOCALDIR="${DIR}/${FAMILY_ID}/${INDIVIDUAL_ID}"

	[ ! -d "$LOCALDIR" ] && $(mkdir -p "$LOCALDIR")

	REL_FILE="${LOCALDIR}/relationship.txt"
	[ ! -f $REL_FILE ] && echo $REL > $REL_FILE

	if [[ $(echo "$FILEPATH" | grep -c "\/mnt" ) == 0 ]]
	then
		echo "Skipping, invalid file path $FILEPATH for $FAMILY_ID $INDIVIDUAL_ID" >> $REPORT_FILE
		continue
	fi
			
	NUMFILES=$(find $FILEPATH -name "*fastq.gz" -or -name "read.bam" | wc -l)
	if [[ $NUMFILES > 1 ]]
	then
		echo "More than one sequence data file in $FILEPATH" >> $REPORT_FILE
		continue
	fi

	for i in $(find $FILEPATH -name "*fastq.gz" -or -name "read.bam")
	do
		TARGETNAME="${LOCALDIR}/${INFO}.bam"

		if [[ ! -f $TARGETNAME ]]
		then
			echo "Linking $i to $TARGETNAME" >> $REPORT_FILE
  			$(ln -s $i $TARGETNAME)
		else
			echo "Already have $TARGETNAME"  >> $REPORT_FILE
		fi
	done
done

echo "Report can be found in $REPORT_FILE"
