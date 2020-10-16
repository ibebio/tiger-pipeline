#!/bin/bash
#####################################################################
# User settings
#
# Number of local cores on the machine where Snakemake is run
# Not too many are needed, since most jobs are run on the cluster
LOCAL_CORES=10

# Allowed number of parallel cluster jobs
JOBS=50

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
snakemake \
    --local-cores ${LOCAL_CORES} \
    --jobs ${JOBS} \
    --use-conda \
    --profile sge \
    ${SNAKEMAKE_ARGS}
