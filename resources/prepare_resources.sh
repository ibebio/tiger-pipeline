#!/bin/bash

# Unzip the reference and index files, if it is the first time the pipeline is run
# Call from the root pipeline directory with bash resources/prepare_resources.sh

cd resources
shopt -s nullglob
FILES=( *.fa*.gz )
if (( ${#FILES[@]} )); then
  for FILE in "${FILES[@]}" ; do
    echo "unzipping ${FILE}"
    gunzip "${FILE}"
  done
fi
shopt -u nullglob
# Create symbolic link for .dict index file
if [[ ! -f Arabidopsis_thaliana.TAIR10.dna.toplevel.dict ]] ; then
   ln -s Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.dict Arabidopsis_thaliana.TAIR10.dna.toplevel.dict
fi
cd ..
