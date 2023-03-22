#!/bin/bash
#####################################################################
# User settings
#
# Number of local cores on the machine where Snakemake is run
# Not too many are needed, since most jobs are run on the cluster
LOCAL_CORES=10

# Allowed number of parallel cluster jobs
JOBS=50



if [[ "$1" != "local" ]] && [[ "$1" != "sge" ]] ; then
    echo "Usage: $0 local|sge [additional Snakemake arguments]"
    echo "The first argument is required, and is either the string 'local' or 'sge'."
    echo "It determines whether the workflow should be run on the local machine (local"
    echo "or on the SGE cluster (sge)."
    echo "Any additional aguments are passed on to Snakemake"
    exit 1
fi

SNAKEMAKE_PROFILE=$1
shift

# Any additional arguments to the script are passed on to Snakemake
# e.g. --forceall would force to rerun all rules and recreate all
#files.
SNAKEMAKE_ARGS="$@"
#
#####################################################################


PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Create environments first, since there is no internet
# connectivity on the cluster nodes
if [[ ! -d ${PIPELINE_DIR}/conda/ ]] ; then
    snakemake --use-conda \
      --conda-create-envs-only \
      --cores 4 \
      --conda-frontend mamba
fi


# Unzip the reference and index files, if it is the first time the pipeline is run
cd resources
for FILE in *.fa*.gz ; do
    echo "unzipping ${FILE}"
    gunzip ${FILE}
done
# Create symbolic link for .dict index file
if [[ ! -f Arabidopsis_thaliana.TAIR10.dna.toplevel.dict ]] ; then
   ln -s Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.dict Arabidopsis_thaliana.TAIR10.dna.toplevel.dict
fi
cd ..

# Run the workflow
if [[ "${SNAKEMAKE_PROFILE}" == "sge" ]] ; then
    snakemake \
	--local-cores ${LOCAL_CORES} \
	--jobs ${JOBS} \
	--use-conda \
	--profile sge \
  --conda-frontend mamba \
	${SNAKEMAKE_ARGS}
else
    snakemake \
	--jobs ${LOCAL_CORES} \
	--use-conda \
  --conda-frontend mamba \
	${SNAKEMAKE_ARGS}
fi
