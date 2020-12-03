#!/bin/bash
# Prepare input files for the TIGER pipeline
TIGER_SCRIPTS_DIR=$1
WORKDIR=$2
TAB_FILE=$3
MARKER_FILE_COMPLETE=$4
MARKER_FILE_CORRECTED=$5
COMPLETE_FILE=$6
CORRECTED_FILE=$7
SAMPLE_NAME=$8

mkdir -p ${WORKDIR}

#get rid of mito and cp reads; create columns for read counts for ref allele and alt allele
# chr pos ref ref_count alt alt_count
awk '{if ($1 ~ /^[1-5]*$/) print $1 "\t" $2 "\t" $5 "\t" $3 "\t" $6 "\t" $4}' ${TAB_FILE} > ${WORKDIR}/${SAMPLE_NAME}.input.temp

#generate "complete" file
perl ${TIGER_SCRIPTS_DIR}/get_subset.pl ${WORKDIR}/${SAMPLE_NAME}.input.temp 1,2 ${MARKER_FILE_COMPLETE} 1,2 0 > ${WORKDIR}/${SAMPLE_NAME}.input.complete.temp

#generate "corrected" file for TIGER)
perl ${TIGER_SCRIPTS_DIR}/get_subset.pl ${WORKDIR}/${SAMPLE_NAME}.input.temp 1,2 ${MARKER_FILE_CORRECTED} 1,2 0 > ${WORKDIR}/${SAMPLE_NAME}.input.corrected.temp

# REMOVE LINES WITH missing ref allele count "." (missing GT), THEN CHANGE REMAINING "." to "0" (no ALT allele call, hom REF), remove lines with bases-entries in count coluimns, restrict to SNPs (allele length =1)
cat ${WORKDIR}/${SAMPLE_NAME}.input.complete.temp | awk '{ if (($4+$6 != "0") && ($4 ~ /^[0-9]*$/) && (length($3)==1) && (length($5)==1)) { print $0}}' | sed 's/\./0/g' > ${COMPLETE_FILE}
cat ${WORKDIR}/${SAMPLE_NAME}.input.corrected.temp | awk '{ if (($4+$6 != "0") && ($4 ~ /^[0-9]*$/) && (length($3)==1) && (length($5)==1)) { print $0}}' | sed 's/\./0/g' > ${CORRECTED_FILE}
