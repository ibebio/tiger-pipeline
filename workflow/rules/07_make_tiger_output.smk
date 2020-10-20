localrules: tiger_create_genotype_plots, tiger_create_path_file, tiger_create_rqtl_input


rule tiger_create_genotype_plot:
    input:
        rough_co_breaks="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.breaks.txt",
        rough_co_breaks_refined="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.refined.breaks.txt",
        breaks_refined_corrected="results/tiger_analysis/F2.{crossing_id}/rough_co_breaks_refined_corrected/{f2_sample}.corrected.refined.breaks.txt",
        frequencies="results/tiger_analysis/F2.{crossing_id}/allele_frequencies/{f2_sample}.input.corrected.frequencies_bmm.txt",
        sliding_window_breaks="results/tiger_analysis/F2.{crossing_id}/hmm_model/{f2_sample}_sliding_window.breaks.txt"
    output:
        plot="results/plots/F2.{crossing_id}/{f2_sample}.pdf",
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        sample=lambda wildcards: wildcards.f2_sample,
        rough_co=lambda wildcards: "results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.txt".format(crossing_id=wildcards.crossing_id, f2_sample=wildcards.f2_sample)
    resources:
    log:
        "results/logs/tiger_create_genotype_plot/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/r.yaml"
    shell:
        """
        # TODO Check whats the problem was with this
        R --slave --vanilla --args \
        {params.sample} \
        {output.plot} \
        {params.rough_co} \
        {input.rough_co_breaks} \
        {input.breaks_refined_corrected} \
        {input.frequencies} \
        {input.sliding_window_breaks} < {params.tiger_scripts_dir}/plot_genotyping.R
        """

rule tiger_create_path_file:
    input:
        corrected_refined_breaks_files=get_corrected_refined_breaks_files
    output:
        breaks_summary="results/tiger_analysis/F2.{crossing_id}/{crossing_id}.breaks_paths.txt",
    params:
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_create_path_file/{crossing_id}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        if [[ -f {output.breaks_summary} ]] ; then \
          rm {output.breaks_summary} ; \
        fi ; \
        for FILE in {input.corrected_refined_breaks_files} ; do \
          SAMPLE=$(basename $FILE |sed 's/.corrected.refined.breaks.txt$//') ; \
          echo "$SAMPLE\t$FILE" >> {output.breaks_summary} ; \
        done
        """


rule tiger_create_rqtl_input:
    input:
        breaks_summary="results/tiger_analysis/F2.{crossing_id}/{crossing_id}.breaks_paths.txt",
    output:
        rqtl_csv="results/tables/rqtl/{crossing_id}.rqtl.csv"
    params:
        tiger_chrom_size_file=config["tiger"]["chromosome_size_file"]
    log:
        "results/logs/tiger_create_rqtl_input/{crossing_id}.log"
    shell:
        """
        env python workflow/scripts/tiger_to_rqtl.py \
           --samples_file {input.breaks_summary} \
           --chrom_size_file {params.tiger_chrom_size_file} \
           --output {output.rqtl_csv} 2> {log}
        """
    
rule tiger_create_genotype_plots_all:
    input:
        plots=get_plot_files_for_crossing_id
    output:
        plot_marker=touch("results/plots/{crossing_id}.done")
