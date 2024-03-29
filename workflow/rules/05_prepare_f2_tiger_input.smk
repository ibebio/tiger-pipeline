# Rule prepare_corrected_inputs_f2:
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

localrules: prepare_tiger_inputs_f2, prepare_tiger_inputs_f2_all,qc_pipeline_overview_image

rule prepare_tiger_inputs_f2:
    input:
        tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tsv",
	    marker_file_complete="results/markers/complete/{crossing_id}.SNP.biallelic.complete.txt",
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

# QC pipeline visualization rule
rule qc_pipeline_overview_image:
    input:
        flag_f2_inputs="results/tiger_analysis/prepare_tiger_inputs_f2.done",
        flag_rqtl_plots="results/plots/{crossing_id}.done",
        src_parent="results/variants/parental/{crossing_id}.src_parent.txt",
        ref_parent="results/variants/parental/{crossing_id}.ref_parent.txt",
        no_te_no_tlr_vcf="results/variants/parental/isec_parent_lines_noTEnoTLR/{crossing_id}.vcf",
        isec_output_dir=directory("results/variants/parental/isec_parent_lines/isec.{crossing_id}")
    output:
        overview_image="results/qc/TIGER_input_summary_for_crossing_id_{crossing_id}_mqc.png"
    params:
        crossing_id="{crossing_id}"
    threads: 1
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/qc_pipeline_overview.py"
