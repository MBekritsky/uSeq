#!/bin/bash

DIR=$1
LB_SCHEME=$2
GENOME="/mnt/wigclust4/home/bekritsk/genomes/hg19/Mod/hg19.mod.fa"
USEQ_CONFIG="/mnt/wigclust4//home/bekritsk/tools/microsatellite/config.txt"
COUNT=0

[ -z $LB_SCHEME ] && LB_SCHEME="2"

REPORT_FILE="${DIR}/run_$(date +%d%m%Y).txt"
[ -f $REPORT_FILE ] && rm $REPORT_FILE 

for fam in $(find $DIR -maxdepth 1 -mindepth 1 -type d -name "auSSC*")
do
	FAMID=$(basename $fam)
	for person in $(find $fam -maxdepth 1 -name "SSC*" -type d)
	do
		PERSONID=$(basename $person)
		PROJECT_ID="simons/${FAMID}/${PERSONID}"
		SGEERR="${person}/${PERSONID}.err"
		SGEOUT="${person}/${PERSONID}.out"		
		echo "Running MicroSeq on $PERSONID with project ID $PROJECT_ID" >> $REPORT_FILE
		OUTPUT=$(qsub -V -P Make -q make@wigclust10 -S /bin/bash -v USEQ_CONFIG=$USEQ_CONFIG -e $SGEERR -o $SGEOUT /mnt/wigclust4/home/bekritsk/tools/microsatellite/MicroSeq2.pl --dir $person --genome $GENOME --sge --project_id $PROJECT_ID --lb_scheme $LB_SCHEME)
		echo -e "\t$OUTPUT" >> $REPORT_FILE
		COUNT=$((COUNT + 1))
	done
done

echo "$COUNT individuals submitted"
echo "Report can be found in $REPORT_FILE"
