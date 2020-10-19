localrules: tiger_create_path_file

rule tiger_create_path_file:
    input:
        corrected_refined_breaks_files=get_corrected_refined_breaks_files # TODO create this fct
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
          echo $FILE >> {output.breaks_summary} ; \
        done
        """
	
    
# rule tiger_plot_output:
#     input:
#         corrected_refined_breaks_files=get_corrected_refined_breaks_files # TODO create this fct
#     output:
#         breaks_summary="results/tiger_analysis/F2.{crossing_id}/{crossing_id}.breaks_paths.txt",
#     params:
#     resources:
#         n=1,
#         time=lambda wildcards, attempt: 12 * 59 * attempt,
#         mem_gb_pt=lambda wildcards, attempt: 12 * attempt
#     log:
#         "results/logs/tiger_create_path_file/{crossing_id}.log"
#     conda:
#         "../envs/r.yaml"
#     shell:
#         """
#         if [[ -f {output.breaks_summary} ]] ; then \
#           rm {output.breaks_summary} ; \
#         fi \ ;
#         for $FILE in {input.corrected_refined_breaks_files} ; do \
#           echo ${FILE} >> {output.breaks_summary} ; \
#         done
#         """


# R --slave --vanilla --args \
#     TEST \
#     TEST.pdf \
#     F2-16-1-004.rough_COs.breaks.txt \
#     F2-16-1-004.rough_COs.refined.breaks.txt \
#     F2-16-1-004.corrected.refined.breaks.txt \
#     F2-16-1-004.input.corrected.frequencies_bmm.txt \
#     F2-16-1-004_sliding_window.breaks.txt < /ebio/abt6_projects9/At_dFLC/code/TIGER_scripts/plot_genotyping.R


