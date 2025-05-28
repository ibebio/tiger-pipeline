# Snakemake workflow: TIGER analysis pipeline



[![Snakemake](https://img.shields.io/badge/snakemake-≥5.21.0-brightgreen.svg)](https://snakemake.bitbucket.io)


## Authors

* Ulrich Lutz
* Snakemake workflow: Ilja Bezrukov

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors of this pipeline and the authors of the 'GBS with sparse coverage using Trained Individual
GenomE Reconstruction (TIGER)' workflow: https://www.g3journal.org/content/5/3/385.short
### Step 1: Obtain a copy of this workflow
Clone the repository into the place where you want to perform the data analysis.
For convenience, and reproducibility, the original TIGER scripts with small bugfixes from the paper of Rowan et al are included.

For reference, the original forked repository of the TIGER scripts with the fixes can be found at (https://github.com/ibebio/TIGER_Scripts-for-distribution)

### Step 2: Configure workflow
Change to the freshly cloned pipeline directory
       
       cd tiger-pipeline

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.csv` to specify your sample setup. The pipeline is shipped with the Arababidopsis thaliana TAIR10 reference.

Run the following command  to make the required scripts executable:
```
$ chmod u+x workflow/scripts/*.*
```

### Step 3: Install Snakemake
Install Snakemake using [mamba](https://github.com/mamba-org/mamba):

	mamba create -c conda-forge -c bioconda -n snakemake snakemake">="5.21.0,"<="6.10.0 python">="3.7

Mamba is an alternative package manager for the conda ecosystem with a much
more reliable dependency resolution and better speed. The pipeline was not tested with Snakemake versions newer than 6.10.0. While the workflow should work, the profile configuration in the following section will need to be adjusted.

In case you have conda available, but not mamba, run 

  conda install -n base -c conda-forge mamba

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).



### Step 4: Execute workflow

##### Basic workflow for any compute environment

If using the provided Arabidopsis thaliana reference, run
```
bash resources/prepare_resources.sh
```
in the pipeline directory before running the pipeline.

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores. 

Then, run the workflow via

    snakemake --use-conda --jobs 100

The number of jobs can be adjusted as required. Additional arguments
for Snakemake to use cluster scheduling can be supplied, please check the Snakemake documentation on how to do this.

To run it in a cluster environment, you might need to create all required conda
environments first via

    snakemake --use-conda --conda-create-envs-only --cores 4


#### Pre-configured workflow for SGE clusters
For the Max-Planck-Institute for Biology Tübingen, set up your SGE cluster profile as follows:

```
git clone https://github.com/ibebio/snakemake_profiles.git
cd snakemake_profiles
mkdir -p ~/.config/snakemake/
chmod u+x sge/*.py
cp -r sge ~/.config/snakemake/
```
    
Activate the conda environment:

    mamba activate snakemake

Test your configuration by performing a dry-run via

    snakemake -n


A helper script `run-workflow.sh`, is included to conveniently run the
workflow, either locally or on the cluster:

	./run-workflow.sh sge

would run the pipeline on the SGE cluster, as set up previously.

	./run-workflow.sh local

would run it on the local maschine.

To customize how many cores and jobs are used, you can either modify
the `run-workflow.sh` script or run the commands required to run the
workflow by hand, as described below.

To clean up all output files and conda environments to rerun the workflow from
scratch, the helper script `clean-all.sh` is included.


### Step 5: Investigate results

All output is stored in the `results/` subfolder.
Logs for each step are stored in `logs/`.

The `workflow/` folder contains the Snakemake files and scripts that are needed to run the workflow.
It does not need to be changed unless the workflow has to be modifed.

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with bugfixes or new developments from the upstream repository, do the following:

1. At the very least, your config files will be different, compared to the example ones from upstream. Therefore, they need to be secured before obtaining the upstream copy: `git stash`
2. Obtain the updates from the Github repository: `git pull`
3. Restore your modifications to the config files: `gut stash pop`

The above steps assume that you did not modify any parts of the workflow, except the config files. If the config format has changed, you might need to update them.