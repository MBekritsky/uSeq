#!/bin/bash

DATA_PATH="/mnt/wigclust4/home/bekritsk/tools/microsatellite"

FILE_ID=$1
TARGET_DIR=$2
SUFFIX=$3
ALN_GENOME=$4
FILE_PATH=$5

FASTQ_FILE_1="${FILE_PATH}s_${FILE_ID}_1_sequence.${SUFFIX}"
FASTQ_FILE_2="${FILE_PATH}s_${FILE_ID}_2_sequence.${SUFFIX}"

echo "Processing ${FASTQ_FILE_1} and ${FASTQ_FILE_2}"

SAI_FILE_1="${TARGET_DIR}s_${FILE_ID}_1_sequence.sai"
SAI_FILE_2="${TARGET_DIR}s_${FILE_ID}_2_sequence.sai"

/mnt/wigclust1/data/software/bwa/default/bwa aln -t 4 "$ALN_GENOME" "$FASTQ_FILE_1" > "$SAI_FILE_1" &
PID_1=$!
/mnt/wigclust1/data/software/bwa/default/bwa aln -t 4 "$ALN_GENOME" "$FASTQ_FILE_2" > "$SAI_FILE_2" &
PID_2=$!

echo $PID_1
echo $PID_2

wait $PID_1
wait $PID_2

echo "Completed bwa aln on $FASTQ_FILE_1 and $FASTQ_FILE_2 with reference genome $ALN_GENOME"

SAMPE_OUT_FILE="${TARGET_DIR}s_${FILE_ID}_sequence.sam"

echo "Generating sam file from $SAI_FILE_1 and $SAI_FILE_2 using sampe"
/data/software/bwa/default/bwa sampe "$ALN_GENOME" "$SAI_FILE_1" "$SAI_FILE_2" "$FASTQ_FILE_1" "$FASTQ_FILE_2" > "$SAMPE_OUT_FILE"
echo "Completed bwa sampe on $SAI_FILE_1 and $SAI_FILE_2"

$(rm $SAI_FILE_1)
$(rm $SAI_FILE_2)

BAM_FILE="${TARGET_DIR}s_${FILE_ID}_sequence.bam"

/data/software/bwa/default/bwa view -bS "$SAMPE_OUT_FILE" -o "$BAM_FILE"
