#!/bin/bash

# user input
if [ "$#" -lt 1 ]; then
    echo "Usage: clean-all.sh [preview|delete]"
    echo "Description: delete all snakemake-generated files"
    echo " preview => preview files to delete"
    echo " delete => delete files"
    echo "NOTE: this only deletes files that snakemake knows about, but that's all that's needed to fully restart the snakemake pipeline"
    exit
fi

PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

FILES=`snakemake --summary --rerun-incomplete | tail -n+2 | cut -f1`
if [ $1 == "preview" ]; then
    echo "#-- Files to delete --#"
    printf '%s\n' "${FILES[@]}"
    ls ${PIPELINE_DIR}/.snakemake/conda* 2> /dev/null
elif [ $1 == "delete" ]; then
    echo "#-- Deleting the following files --#"
    printf '%s\n' "${FILES[@]}"
    rm -rf $FILES
    if [ -d ${PIPELINE_DIR}/.snakemake/conda ] ; then
      echo "# -- Deleting conda environments --#"
      rm -rf ${PIPELINE_DIR}/.snakemake/conda* 
    fi
else
    echo "$1 not recoginized"
fi 
