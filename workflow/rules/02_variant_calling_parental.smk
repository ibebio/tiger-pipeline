localrules: estimate_filtering_params_parental_corrected
rule call_variants_parental:
    input:
        bam="results/rmdup/{parental_sample}.rmdup.bam"
    output:
        vcf="results/variants/parental/raw/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["gatk_options"]["java_options"],
    threads: 4
    resources:
        n=4,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/call_variants_parental/{parental_sample}.log"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "{params.java_options}" HaplotypeCaller \
            -R {params.index} \
            -I {input.bam} \
            -O {output.vcf} 2> {log}
        """


rule filter_variants_parental_complete:
    input:
        vcf="results/variants/parental/raw/{parental_sample}.vcf",
    output:
        snps=temp("results/variants/parental/raw/{parental_sample}.snps.vcf"),
        indels=temp("results/variants/parental/raw/{parental_sample}.indels.vcf"),
        filtered_snps=temp("results/variants/parental/filtered_complete/{parental_sample}.snps.vcf"),
        filtered_indels=temp("results/variants/parental/filtered_complete/{parental_sample}.indels.vcf"),
        filtered_vcf="results/variants/parental/filtered_complete/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        snp_filter=config["variant_filtering_parental"]["snps_complete_filter"],
        indel_filter=config["variant_filtering_parental"]["indels_complete_filter"],
        java_options=config["gatk_options"]["java_options"],
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/filter_variants_parental_complete/{parental_sample}.log"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "{params.java_options}" SelectVariants \
             -R {params.index} \
             -V {input.vcf} \
             --select-type-to-include SNP \
             -O {output.snps} \
        ; \
        gatk --java-options "{params.java_options}" VariantFiltration \
             -R {params.index} \
             -V {output.snps} \
             --filter-name "snps-hard-filter" \
             --filter-expression "{params.snp_filter}" \
             -O {output.filtered_snps} \
        ; \
        gatk --java-options "{params.java_options}" SelectVariants \
             -R {params.index} \
             -V {input.vcf} \
             --select-type-to-include INDEL \
             -O {output.indels} \
        ; \
        gatk --java-options "{params.java_options}" VariantFiltration \
             -R {params.index} \
             -V {output.indels} \
             --filter-name "indels-hard-filter" \
             --filter-expression "{params.indel_filter}" \
             -O {output.filtered_indels} \
        ; \
        picard MergeVcfs \
               INPUT={output.filtered_snps} \
               INPUT={output.filtered_indels} \
               OUTPUT={output.filtered_vcf} 2> {log}
        """


rule estimate_filtering_params_parental_corrected:
    input:
        vcf=get_parental_raw_variants_all # TODO make the filtering per crossing
    output:
        filter_file="results/variants/parental/snps_indels_corrected_filter.txt",
        plot="results/plots/parental_filtering.png"
    params:
        snp_indel_filter=config["variant_filtering_parental"]["snps_indels_corrected_filter"],
        plot_prefix="results/plots/parental_filtering_tmp"
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/estimate_filtering_params_parental_corrected.log"
    conda:
        "../envs/estimate_filtering_params.yaml"
    shell:
        """
        # Prepare input args
        INPUTARGS=$(echo " {input.vcf}" | sed 's/ / --vcf /g')
        # Estimate filtering parameters
        env python workflow/scripts/estimate_parental_filtering_params.py \
            $INPUTARGS \
            --filter_file {output.filter_file} \
            --plot {params.plot_prefix} 2> {log}

        if [[ "{params.snp_indel_filter}" != "auto" ]] ; then
          # Use hard-coded filtering, overwrite filter_file
          echo "{params.snp_indel_filter}" > {output.filter_file}
        fi

        montage {params.plot_prefix}*.png -geometry +0+0 -title "$(cat {output.filter_file})" {output.plot} 2>> {log}
        rm {params.plot_prefix}*.png
        """


rule filter_variants_parental_corrected:
    input:
        vcf="results/variants/parental/filtered_complete/{parental_sample}.vcf",
        filter_file="results/variants/parental/snps_indels_corrected_filter.txt",
    output:
        filtered_vcf="results/variants/parental/filtered_corrected/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["gatk_options"]["java_options"]
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/filter_variants_parental_corrected/{parental_sample}.log"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "{params.java_options}" VariantFiltration \
             -R {params.index} \
             -V {input.vcf} \
             --filter-name "snps-indels-complete-filter" \
             --filter-expression "$(cat {input.filter_file})" \
             -O {output.filtered_vcf} 2> {log}
        """
