# rule prepare_corrected_inputs_f2:
#     input:
# 	tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tabbed.txt",
# 	marker_file_complete=lambda wildcards: "results/markers/complete/{parental_sample}.SNP.biallelic.complete.txt".format(parental_sample=config[wildcards.crossing_id]["parent_ref"]),
# 	marker_file_corrected="results/markers/corrected/{crossing_id}.SNP.private.txt"
#     output:
# 	tig_in_complete="results/tiger_analysis/F2.{crossing_id}/input.complete/{f2_sample}.input.complete",
# 	tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected",
# 	work_dir=temp(directory("results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.workdir"))
#     params:
# 	tiger_scripts_dir=config["tiger"]["scripts_dir"],
# 	sample_name=get_f2_sample_name
#     resources:
#         n=1,
#         time=lambda wildcards, attempt: 12 * 59 * attempt,
#         mem_gb_pt=lambda wildcards, attempt: 12 * attempt
#     log:
#         "results/logs/prepare_corrected_inputs_f2/{f2_sample}.log"
#     conda:
#         "../envs/global.yaml"
#     shell:
#         """
# 	../scripts/prepare_tiger_inputs_f2.sh \
# 	  {params.tiger_scripts_dir} \
# 	  {output.work_dir} \
#           {input.tab} \
#           {input.marker_file_complete} \
#           {input.marker_file_corrected} \
#           {output.tig_in_complete} \
#           {output.tig_in_corrected} \
#           {params.sample_name}
#         """

localrules: prepare_tiger_inputs_f2, prepare_tiger_inputs_f2_all

rule prepare_tiger_inputs_f2:
    input:
        tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tsv",
	marker_file_complete=lambda wildcards: "results/markers/complete/{parental_sample}.SNP.biallelic.complete.txt".format(parental_sample=[c["parent_a"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
	marker_file_corrected="results/markers/corrected/{crossing_id}.SNP.private.txt"
    output:
        tig_in_complete="results/tiger_analysis/F2.{crossing_id}/input.complete/{f2_sample}.input.complete",
	    tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected",
	    work_dir=temp(directory("results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.workdir"))
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
	sample_name=lambda wildcards : wildcards.f2_sample
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/prepare_corrected_inputs_f2/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/global.yaml"
    shell:
        """
	workflow/scripts/prepare_tiger_inputs_f2.sh \
	  {params.tiger_scripts_dir} \
	  {output.work_dir} \
          {input.tab} \
          {input.marker_file_complete} \
          {input.marker_file_corrected} \
          {output.tig_in_complete} \
          {output.tig_in_corrected} \
          {params.sample_name}
        """



# Aggregate rule
rule prepare_tiger_inputs_f2_all:
    input:
        # get_f2_biallelic_snps_corrected_all,
        get_f2_tiger_inputs_all
    output:
        # flag_vcfs=touch("results/variants/f2/monomorphic/call_variants_f2.done"),
        flag_inputs=touch("results/tiger_analysis/prepare_tiger_inputs_f2.done"),
