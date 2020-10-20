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
    snakemake --use-conda --conda-create-envs-only --cores 4
fi

# Run the workflow

if [[ "${SNAKEMAKE_PROFILE}" == "sge" ]] ; then
    snakemake \
	--local-cores ${LOCAL_CORES} \
	--jobs ${JOBS} \
	--use-conda \
	--profile sge \
	${SNAKEMAKE_ARGS}
else
    snakemake \
	--jobs ${LOCAL_CORES} \
	--use-conda \
	${SNAKEMAKE_ARGS}
fi

