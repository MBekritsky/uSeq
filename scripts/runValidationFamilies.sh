#!/bin/bash

DIR=$1
STAGE=$2
LB_SCHEME=$3
GENOME="/mnt/wigclust4/home/bekritsk/genomes/hg19/Mod/hg19.mod.fa"
COUNT=0

[ -z $LB_SCHEME ] && LB_SCHEME="2"
[ -z $STAGE ] && STAGE="scan"

REPORT_FILE="${DIR}/run_$(date +%d%m%Y).txt"
[ -f $REPORT_FILE ] && rm $REPORT_FILE 

for fam in $(find $DIR -maxdepth 1 -mindepth 1 -type d)
do
	FAMID=$(basename $fam)
	for person in $(find $fam -maxdepth 1 -mindepth 1 -type d)
	do
		PERSONID=$(basename $person)
		PROJECT_ID="simons/validation/${FAMID}/${PERSONID}"
		SGEERR="${person}/${PERSONID}.err"
		SGEOUT="${person}/${PERSONID}.out"		
		echo "Running MicroSeq on $PERSONID with project ID $PROJECT_ID" >> $REPORT_FILE
		OUTPUT=$(qsub -V -P Make -q make@wigclust10 -S /bin/bash -b y -e $SGEERR -o $SGEOUT /mnt/wigclust4/home/bekritsk/tools/microsatellite/MicroSeq2.pl \
		--dir $person --genome $GENOME --sge --project_id $PROJECT_ID --lb_scheme $LB_SCHEME --keep_duplicates --stage $STAGE)
		echo -e "\t$OUTPUT" >> $REPORT_FILE
		COUNT=$((COUNT + 1))
	done
done

echo "$COUNT individuals submitted"
echo "Report can be found in $REPORT_FILE"
