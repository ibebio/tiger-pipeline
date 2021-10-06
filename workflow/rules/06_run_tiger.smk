localrules: tiger_basecaller, tiger_allele_freq_estimator, tiger_beta_mixture_model, tiger_prepare_hmm, tiger_calc_transmission_emission_prob, tiger_run_hmm, tiger_estimate_recombination_breakpoints, tiger_refine_recombination_breakpoints, tiger_smooth_out_breaks

rule tiger_basecaller:
    input:
        tig_input_computation_done="results/tiger_analysis/prepare_tiger_inputs_f2.done", # This is not strictly needed, but added this to force the pipeline to create all input files first. Otherwise, it will only run with the 2 parent samples and one F2 samples until the 1st chekpoint
        tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected",
    output:
        basecall="results/tiger_analysis/F2.{crossing_id}/basecalling/{f2_sample}.input.corrected.basecall.txt"
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt,
        java_options=lambda wildcards, attempt: "-Xmx{}g -Xms{}g".format(8 * attempt, 8 * attempt)
    log:
        "results/logs/tiger_basecaller/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        java {resources.java_options} -jar {params.tiger_scripts_dir}/base_caller.jar \
          -r {input.tig_in_corrected} \
          -o {output.basecall} \
          -n bi 2> {log}
        """


rule tiger_allele_freq_estimator:
    input:
        basecall="results/tiger_analysis/F2.{crossing_id}/basecalling/{f2_sample}.input.corrected.basecall.txt",
        tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected"
    output:
        frequencies="results/tiger_analysis/F2.{crossing_id}/allele_frequencies/{f2_sample}.input.corrected.frequencies_bmm.txt"
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt,
        java_options=lambda wildcards, attempt: "-Xmx{}g -Xms{}g".format(8 * attempt,8 * attempt)
    log:
        "results/logs/tiger_allele_freq_estimator/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        java {resources.java_options} -jar {params.tiger_scripts_dir}/allele_freq_estimator.jar \
          -r {input.tig_in_corrected} \
          -o {output.frequencies} \
          -n bi -w 1000 2> {log}
        """


checkpoint tiger_beta_mixture_model:
    input:
        frequencies="results/tiger_analysis/F2.{crossing_id}/allele_frequencies/{f2_sample}.input.corrected.frequencies_bmm.txt",
        tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected"
    output:
        bmm_run_done="results/tiger_analysis/F2.{crossing_id}/beta_mixture_models/{f2_sample}.bmm.intersections.txt.done"
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        minimal_corrected_input_number=config["tiger"]['minimal_corrected_input_number'],
        bmm="results/tiger_analysis/F2.{crossing_id}/beta_mixture_models/{f2_sample}.bmm.intersections.txt"
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_beta_mixture_model/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/r.yaml"
    shell:
        """
        MARKER_COUNT=$(cat {input.tig_in_corrected} | wc -l)
        if [[ $MARKER_COUNT -gt {params.minimal_corrected_input_number} ]] ; then
          set +e
          Rscript --vanilla {params.tiger_scripts_dir}/beta_mixture_model.R \
            {input.frequencies} \
            {params.bmm}
          exitcode=$?
          echo $exitcode > {output.bmm_run_done}
        else
          echo 1 > {output.bmm_run_done}
          echo "Skipping creating {params.bmm} due to low marker count: $MARKER_COUNT"
        fi
        """

        
# 5 Prepare files for HMM probability estimation using the BASECALLER output and the output of the beta mixture model
rule tiger_prepare_hmm:
    input:
        basecall="results/tiger_analysis/F2.{crossing_id}/basecalling/{f2_sample}.input.corrected.basecall.txt",
        tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected",      
        frequencies="results/tiger_analysis/F2.{crossing_id}/allele_frequencies/{f2_sample}.input.corrected.frequencies_bmm.txt",
        bmm="results/tiger_analysis/F2.{crossing_id}/beta_mixture_models/{f2_sample}.bmm.intersections.txt"
    output:
        hmm_prob="results/tiger_analysis/F2.{crossing_id}/hmm_probabilities/{f2_sample}.hmm_probabilities.txt",
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        tiger_chrom_size_file=config["tiger"]["chromosome_size_file"],
        sample=lambda wildcards: wildcards.f2_sample
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_prepare_hmm/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/r.yaml"
    shell:
        """
        perl {params.tiger_scripts_dir}/prep_prob.pl \
          -s {params.sample} \
          -m {input.tig_in_corrected} \
          -b {input.basecall} \
          -c {params.tiger_chrom_size_file} \
          -o {output.hmm_prob} 2> {log}
        """

# 6 Calculate transmission and emission probabilities for the HMM
rule tiger_calc_transmission_emission_prob:
    input:
        basecall="results/tiger_analysis/F2.{crossing_id}/basecalling/{f2_sample}.input.corrected.basecall.txt",
        frequencies="results/tiger_analysis/F2.{crossing_id}/allele_frequencies/{f2_sample}.input.corrected.frequencies_bmm.txt",
        bmm="results/tiger_analysis/F2.{crossing_id}/beta_mixture_models/{f2_sample}.bmm.intersections.txt",
        hmm_prob="results/tiger_analysis/F2.{crossing_id}/hmm_probabilities/{f2_sample}.hmm_probabilities.txt"
    output:
        hmm_model="results/tiger_analysis/F2.{crossing_id}/hmm_model/{f2_sample}_hmm_model",
        sliding_window_breaks="results/tiger_analysis/F2.{crossing_id}/hmm_model/{f2_sample}_sliding_window.breaks.txt"
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        tiger_chrom_size_file=config["tiger"]["chromosome_size_file"],
        prefix=lambda wildcards: "results/tiger_analysis/F2.{crossing_id}/hmm_model/{f2_sample}".format(crossing_id=wildcards.crossing_id, f2_sample=wildcards.f2_sample)
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_calc_tr_em_prob/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        perl {params.tiger_scripts_dir}/hmm_prob.pl \
          -s {input.frequencies} \
          -p {input.hmm_prob} \
          -o {params.prefix} \
          -a {input.bmm} \
          -c {params.tiger_chrom_size_file} 2> {log}
        """


# 7 Run the HMM
rule tiger_run_hmm:
    input:
        basecall="results/tiger_analysis/F2.{crossing_id}/basecalling/{f2_sample}.input.corrected.basecall.txt",
        hmm_prob="results/tiger_analysis/F2.{crossing_id}/hmm_probabilities/{f2_sample}.hmm_probabilities.txt",
        hmm_model="results/tiger_analysis/F2.{crossing_id}/hmm_model/{f2_sample}_hmm_model"
    output:
        hmm_output="results/tiger_analysis/F2.{crossing_id}/hmm_output/{f2_sample}.hmm.out.txt",
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],        
        sample=lambda wildcards: wildcards.f2_sample
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt,
        java_options=lambda wildcards, attempt: "-Xmx{}g -Xms{}g".format(8 * attempt,8 * attempt)
    log:
        "results/logs/tiger_run_hmm/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        java {resources.java_options} -jar {params.tiger_scripts_dir}/hmm_play.jar \
          -r {input.basecall} \
          -o {output.hmm_output} \
          -t bi \
          -z {input.hmm_model} 2> {log}
        """

# 8 Get rough estimate of recombination breakpoint positions
rule tiger_estimate_recombination_breakpoints:
    input:
        hmm_output="results/tiger_analysis/F2.{crossing_id}/hmm_output/{f2_sample}.hmm.out.txt",
        tig_in_corrected="results/tiger_analysis/F2.{crossing_id}/input.corrected/{f2_sample}.input.corrected",
    output:
        #rough_co="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.txt",
        rough_co_breaks="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.breaks.txt"
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        tiger_chrom_size_file=config["tiger"]["chromosome_size_file"],
        sample=lambda wildcards: wildcards.f2_sample,
        rough_co=lambda wildcards: "results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.txt".format(crossing_id=wildcards.crossing_id, f2_sample=wildcards.f2_sample)
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_estimate_recombination_breakpoints/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        perl {params.tiger_scripts_dir}/prepare_break.pl \
          -s {params.sample} \
          -m {input.tig_in_corrected} \
          -b {input.hmm_output} \
          -c {params.tiger_chrom_size_file} \
          -o {params.rough_co}
        """

        
# 9 Refine recombination breaks
rule tiger_refine_recombination_breakpoints:
    input:
        rough_co_breaks="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.breaks.txt",
        tig_in_complete="results/tiger_analysis/F2.{crossing_id}/input.complete/{f2_sample}.input.complete",
    output:
        rough_co_breaks_refined="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.refined.breaks.txt",
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
        tiger_chrom_size_file=config["tiger"]["chromosome_size_file"],
        sample=lambda wildcards: wildcards.f2_sample
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_refine_recombination_breakpoints/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        perl {params.tiger_scripts_dir}/refine_recombination_break.pl \
	  {input.tig_in_complete} \
	  {input.rough_co_breaks}
        """

# 10 Smooth out breaks
# output .refined.breaks.txt
rule tiger_smooth_out_breaks:
    input:
        rough_co_breaks_refined="results/tiger_analysis/F2.{crossing_id}/rough_co/{f2_sample}.rough_COs.refined.breaks.txt",
    output:
        breaks_refined_corrected="results/tiger_analysis/F2.{crossing_id}/rough_co_breaks_refined_corrected/{f2_sample}.corrected.refined.breaks.txt",
    params:
        tiger_scripts_dir=config["tiger"]["scripts_dir"],
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/tiger_smooth_out_breaks/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/tiger.yaml"
    shell:
        """
        perl {params.tiger_scripts_dir}/breaks_smoother.pl \
	      -b {input.rough_co_breaks_refined} \
          -o {output.breaks_refined_corrected}
        """

# use for mapping, final output: .corrected.refined.breaks.txt
