from snakemake.utils import validate
import pandas as pd
import itertools

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
# validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t|;|,").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
# validate(samples, schema="../schemas/samples.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index),
    # TODO FIX parental_sample="|".join(list(itertools.chain(*[c["parent_ref"] for c in config["crossings"]]))) + "|" + "|".join(list(itertools.chain(*[c["parent_alt"] for c in config["crossings"]]))) ,
    f2_sample="|".join(list(itertools.chain(*[c["f2_samples"] for c in config["crossings"]])))


##### Helper functions #####
def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_all_mapped_files(wildcards):
    """ return the names of all rmdup files (output of 01_mapping.smk) """
    mapped_files = ['results/rmdup/{}.rmdup.bam'.format(s) for s in samples.index]
    return mapped_files


def get_parental_raw_variants_all(wildcards):
    """ return the names of all raw parental variants """
    mapped_files = []
    for c in config["crossings"]:
        mapped_files.append(
            "results/variants/parental/raw/{}.vcf".format(c["parent_ref"]))
        mapped_files.append(
            "results/variants/parental/raw/{}.vcf".format(c["parent_alt"]))
    return mapped_files



def get_parental_biallelic_snps_corrected_all(wildcards):
    """ return the names of all rmdup files (output of 01_mapping.smk) """
    mapped_files = []
    for c in config["crossings"]:
        mapped_files.append(
            "results/variants/parental/biallelic_snps_corrected/{}.vcf".format(c["parent_ref"]))
        mapped_files.append(
            "results/variants/parental/biallelic_snps_corrected/{}.vcf".format(c["parent_alt"]))
    return mapped_files


# def get_f2_biallelic_snps_corrected_all(wildcards):
#     """ return the names of f2 vcfs (output of 04_variant_calling_f2.smk) """
#     mapped_files = []
#     for c in config["crossings"]:
#         mapped_files += \
#             ["results/variants/f2/monomorphic/{crossing_id}{f2_sample}.monomorphic.vcf.gz".format(
#                 crossing_id=wildcards.crossing_id, f2_sample=f) for f in config[c]["f2_samples"]]
#     return mapped_files


def get_f2_tiger_inputs_all(wildcards):
    outputs = []
    for c in config["crossings"]:
        outputs += [
                "results/tiger_analysis/F2.{crossing_id}/input.complete/{f2_sample}.input.complete".format(crossing_id=c["id"], f2_sample=f) for f in c["f2_samples"]]
        outputs += [
                "results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected".format(crossing_id=c["id"], f2_sample=f) for f in c["f2_samples"]]

    return outputs
    

def get_parental_complete_markers_all(wildcards):
    """ return the names of all complete marker files (output of 01_mapping.smk) """
    marker_files = []
    for c in config["crossings"]:
        marker_files.append(
            "results/markers/complete/{sample}.SNP.biallelic.complete.txt".format(sample=c["parent_ref"]))
        marker_files.append(
            "results/markers/complete/{sample}.SNP.biallelic.complete.txt".format(sample=c["parent_alt"]))
    return marker_files


# def get_f2_sample_name(wildcards): # Do we really need this?
#     """Get the sample name for the current f2 sample. Needed for params"""
#     return config["crossings"][wildcards.crossing_id][wildcards.f2_sample]


def get_corrected_refined_breaks_files(wildcards):
    """ Get the filenames of all breaks_files for a crossing id """
    f2_sample_names = [c["f2_samples"] for c in config["crossings"] if c["id"] == wildcards.crossing_id]
    assert len(f2_sample_names) == 1
    f2_sample_names = f2_sample_names[0]
    corrected_refined_breaks_files = []
    for f2_sample in f2_sample_names:
        with checkpoints.tiger_beta_mixture_model.get(f2_sample=f2_sample, crossing_id=wildcards.crossing_id).output[0].open() as f:
            # Was the beta mixture model successfully generated?
            if f.read().strip() == "0":
                corrected_refined_breaks_files.append(
                    "results/tiger_analysis/F2.{crossing_id}/rough_co_breaks_refined_corrected/{f2_sample}.corrected.refined.breaks.txt".format(basedir=workflow.basedir, f2_sample=f2_sample, crossing_id=wildcards.crossing_id)
                    )
    return corrected_refined_breaks_files


def get_plot_files_for_crossing_id(wildcards):
    """ Get the filenames of all plot files for a crossing id """
    f2_sample_names = [c["f2_samples"] for c in config["crossings"] if c["id"] == wildcards.crossing_id]
    assert len(f2_sample_names) == 1
    f2_sample_names = f2_sample_names[0]

    plot_files = []
    for f2_sample in f2_sample_names:
        with checkpoints.tiger_beta_mixture_model.get(f2_sample=f2_sample, crossing_id=wildcards.crossing_id).output[0].open() as f:
            # Was the beta mixture model successfully generated?
            if f.read().strip() == "0":
                plot_files.append(
                    "results/plots/F2.{crossing_id}/{f2_sample}.pdf".format(crossing_id=wildcards.crossing_id, f2_sample=f2_sample)
                )
    return plot_files
