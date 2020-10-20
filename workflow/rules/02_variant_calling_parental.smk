rule call_variants_parental:
    input:
        bam="results/rmdup/{parental_sample}.rmdup.bam"
    output:
        vcf="results/variants/parental/raw/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["variant_calling_parental"]["java_options"],
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
        indel_filter=config["variant_filtering_parental"]["indels_complete_filter"]
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
        gatk SelectVariants \
             -R {params.index} \
             -V {input.vcf} \
             --select-type-to-include SNP \
             -O {output.snps} \
        ; \
        gatk VariantFiltration \
             -R {params.index} \
             -V {output.snps} \
             --filter-name "snps-hard-filter" \
             --filter-expression "{params.snp_filter}" \
             -O {output.filtered_snps} \
        ; \
        gatk SelectVariants \
             -R {params.index} \
             -V {input.vcf} \
             --select-type-to-include INDEL \
             -O {output.indels} \
        ; \
        gatk VariantFiltration \
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


rule filter_variants_parental_corrected:
    input:
        vcf="results/variants/parental/raw/{parental_sample}.vcf",
    output:
        filtered_vcf="results/variants/parental/filtered_corrected/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        snp_indel_filter=config["variant_filtering_parental"]["snps_indels_corrected_filter"]
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
        gatk VariantFiltration \
             -R {params.index} \
             -V {input.vcf} \
             --filter-name "snps-indels-complete-filter" \
             --filter-expression "{params.snp_indel_filter}" \
             -O {output.filtered_vcf} \
        """
