#!/bin/bash

PEOPLE_FILE=$1
LB_SCHEME=$3
GENOME="/mnt/wigclust4/home/bekritsk/genomes/hg19/Mod/hg19.mod.fa"
USEQ_CONFIG="/mnt/wigclust4//home/bekritsk/tools/microsatellite/config.txt"
COUNT=0

[ -z $LB_SCHEME ] && LB_SCHEME="2"

DIR=$(dirname $PEOPLE_FILE)
REPORT_FILE="${DIR}/rerun_$(date +%d%m%Y).txt"
[ -f $REPORT_FILE ] && rm $REPORT_FILE

while read line
do
	person=$line
	STAGE=$(echo "$person" | awk '{print $2}')
	DIR=$(echo "$person" | awk '{print $1}')

	[ -z $STAGE ] && STAGE="scan"

	[ $STAGE == "OK" ] && continue
	
	FAMID=$(echo $DIR | cut -f1 | tr '/' '\n' | grep auSSC)
	PERSONID=$(basename $DIR)
	PROJECT_ID="simons/${FAMID}/${PERSONID}"
	SGEERR="${DIR}/${PERSONID}.err"
	SGEOUT="${DIR}/${PERSONID}.out"

	if [ $STAGE == "scan" ]
	then
	 [ -d "${DIR}/SGEERR" ] && rm -r "${DIR}/SGEERR"
	 [ -d "${DIR}/SGEOUT" ] && rm -r "${DIR}/SGEOUT"
	elif [ $STAGE == "profile" ]
	then
	 [ -e "${DIR}/SGEERR/${PERSONID}.md.err" ] && rm -r "${DIR}/SGEERR/${PERSONID}.md.err"
	 [ -e "${DIR}/SGEOUT/${PERSONID}.md.out" ] && rm -r "${DIR}/SGEOUT/${PERSONID}.md.out"
	fi

	echo "Running MicroSeq on $PERSONID at stage $STAGE with project ID $PROJECT_ID" >> $REPORT_FILE
	OUTPUT=$(qsub -V -P Make -q make@wigclust10 -S /bin/bash -b y -v USEQ_CONFIG=$USEQ_CONFIG -e $SGEERR -o $SGEOUT /mnt/wigclust4/home/bekritsk/tools/microsatellite/MicroSeq2.pl --dir $DIR --genome $GENOME --sge --project_id $PROJECT_ID --stage $STAGE --lb_scheme $LB_SCHEME)
#	echo "qsub -V -P Make -q make@wigclust10 -S /bin/bash -b y -e $SGEERR -o $SGEOUT /mnt/wigclust4/home/bekritsk/tools/microsatellite/MicroSeq2.pl --dir $DIR --genome $GENOME --sge --project_id $PROJECT_ID --stage $STAGE --lb_scheme $LB_SCHEME" >> $REPORT_FILE
	echo -e "\t$OUTPUT" >> $REPORT_FILE
	COUNT=$((COUNT + 1))
done < $PEOPLE_FILE

echo "$COUNT individuals submitted"
echo "Report can be found in $REPORT_FILE"
