#!/bin/bash

# Performs
# - mapping via bwa mem
# - add read group in the process
# - sort via samtools

THREADS=$1
PLATFORM=$2
SAMPLE=$3
INDEX=$4
OUTPUT=$5
FQ1=$6
FQ2=$7
LIBRARY=$8
TEMP_BAM=${OUTPUT%.*}.tmp

bwa mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}\tPU:${PLATFORM}" \
    ${INDEX} \
    ${FQ1} \
    ${FQ2} \
    | \
    samtools sort -@ ${THREADS} -o ${OUTPUT} -T ${TEMP_BAM}
# samtools view -@ 9 -Sb - > $OUTPUT
