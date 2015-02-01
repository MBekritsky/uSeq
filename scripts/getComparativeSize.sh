#!/bin/bash
FAMILY=$1


#for FAMILY in $(find . -name "auSSC*")
#do
	for INDIVIDUAL in $(find $FAMILY  -mindepth 1 -maxdepth 1 -type d)
	do
	FASTQ_SIZE=0
	TOTAL_SIZE=0
		for i in $(ls -al ${INDIVIDUAL} | grep fastq | awk '{print $11}')
		do
			FASTQ_SIZE=$(($FASTQ_SIZE + $(ls -al $i | awk '{print $5}')))
		done
		TOTAL_SIZE=$(($(du "${INDIVIDUAL}" | tail -1 | cut -f1)*1000))
	
	echo $INDIVIDUAL
	FASTQ_SIZE_G=`perl -e 'printf("%0.2fG",$ARGV[0]/1e9)' $FASTQ_SIZE`
	TOTAL_SIZE_G=`perl -e 'printf("%0.2fG",$ARGV[0]/1e9)' $TOTAL_SIZE`
	RATIO=`perl -e 'printf("%0.2f%%",(($ARGV[0]/$ARGV[1])*100))' $TOTAL_SIZE $FASTQ_SIZE`
	echo "Sequencing file size: $FASTQ_SIZE_G"
	echo "MicroSeq files size: $TOTAL_SIZE_G"
	echo $RATIO
	done
#done
