rule collect_allelic_counts_f2:
    input:
        bam="results/rmdup/{f2_sample}.rmdup.bam",
        marker_complete_src_vcf="results/variants/parental/biallelic_snps_complete/src.{crossing_id}.vcf"
    output:
        tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.raw.tsv"
    params:
        index=config["ref"]["genome"],
        java_options=config["gatk_options"]["java_options"]
    resources:
        n=1,
        time=lambda wildcards, attempt: 1 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 60 * attempt
    log:
        "results/logs/collect_allelic_counts_f2/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/gatk4.yaml"
    shell:
        """
        gatk --java-options "{params.java_options}" CollectAllelicCounts \
            -R {params.index} \
            -I {input.bam} \
            -L {input.marker_complete_src_vcf} \
            -O {output.tab} 2> {log}
        """


rule filter_allelic_counts_f2:
    input:
        tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.raw.tsv",
    output:
        tab="results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tsv"
    params:
        sd_factor=config["filtering_f2"]["sd_factor"],
        sample=lambda wildcards : wildcards.f2_sample
    resources:
        n=1,
        time=lambda wildcards, attempt: 1 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 60 * attempt
    log:
        "results/logs/filter_allelic_counts_f2/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/estimate_filtering_params.yaml"
    shell:
        """
        env python workflow/scripts/filter_allelic_counts.py \
           --name {params.sample} \
           --input_tab {input.tab} \
           --output_tab {output.tab} \
           --sd_factor {params.sd_factor}
        2> {log}
        """


###############################################################################
# DEBUG RULES for creating of F2 VCFs
###############################################################################
if config['debug']['call_variants_f2']['create_vcf_for_f2'] == 'yes':
    rule dbg_call_variants_f2:
        input:
            bam="results/rmdup/{f2_sample_vc}.rmdup.bam"
        output:
            gvcf=temp("results/variants/f2/gvcf/{crossing_id}.{f2_sample_vc}.g.vcf.gz"),
        params:
            index=config["ref"]["genome"],
            java_options=config["gatk_options"]["java_options"],
        threads: 4
        resources:
            n=4,
            time=lambda wildcards, attempt: 48 * 59 * attempt,
            mem_gb_pt=lambda wildcards, attempt: 12 * attempt
        log:
            "results/logs/dbg_call_variants_f2/{crossing_id}.{f2_sample_vc}.log"
        conda:
            "../envs/gatk4.yaml"
        shell:
            """
            gatk --java-options "{params.java_options}" HaplotypeCaller \
                -ERC GVCF \
                -R {params.index} \
                -I {input.bam} \
                -O {output.gvcf} 2> {log}
            """


    rule dbg_combine_calls_f2:
        input:
            gvcfs=lambda wildcards: expand("results/variants/f2/gvcf/{{crossing_id}}.{f2_sample_vc}.g.vcf.gz",
                f2_sample_vc=[c["f2_samples"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0])

        output:
            gvcf=temp("results/variants/f2/gvcf/{crossing_id}.g.vcf.gz"),
        params:
            index=config["ref"]["genome"],
            java_options=config["gatk_options"]["java_options"],
        threads: 4
        resources:
            n=4,
            time=lambda wildcards, attempt: 24 * 59 * attempt,
            mem_gb_pt=lambda wildcards, attempt: 12 * attempt
        log:
            "results/logs/dbg_combine_calls_f2/{crossing_id}.log"
        conda:
            "../envs/gatk4.yaml"
        shell:
            """
            gatk --java-options "{params.java_options}" CombineGVCFs \
                -V $(echo "{input.gvcfs}" | sed 's/ / -V /g') \
                -R {params.index} \
                -O {output.gvcf} 2> {log}
            """


    rule dbg_genotype_variants_f2:
        input:
            gvcf="results/variants/f2/gvcf/{crossing_id}.g.vcf.gz",
        output:
            vcf="results/variants/raw/{crossing_id}.vcf.gz",
        params:
            index=config["ref"]["genome"],
            java_options=config["gatk_options"]["java_options"],
        threads: 4
        resources:
            n=4,
            time=lambda wildcards, attempt: 24 * 59 * attempt,
            mem_gb_pt=lambda wildcards, attempt: 12 * attempt
        log:
            "results/logs/dbg_genotype_variants_f2/{crossing_id}.log"
        conda:
            "../envs/gatk4.yaml"
        shell:
            """
            gatk --java-options "{params.java_options}" GenotypeGVCFs \
                -V {input.gvcf} \
                -R {params.index} \
                -O {output.vcf} 2> {log}
            """


    rule dbg_filter_variants_f2:
        input:
            vcf="results/variants/raw/{crossing_id}.vcf.gz",
        output:
            snps=temp("results/variants/f2/raw/{crossing_id}.snps.vcf"),
            indels=temp("results/variants/f2/raw/{crossing_id}.indels.vcf"),
            filtered_snps=temp("results/variants/f2/filtered/{crossing_id}.snps.vcf"),
            filtered_indels=temp("results/variants/f2/filtered/{crossing_id}.indels.vcf"),
            filtered_vcf="results/variants/f2/filtered/{crossing_id}.vcf",
        params:
            index=config["ref"]["genome"],
            snp_filter=config["debug"]["call_variants_f2"]["variant_filtering"]["snps_filter"],
            indel_filter=config["debug"]["call_variants_f2"]["variant_filtering"]["indels_filter"],
            java_options=config["gatk_options"]["java_options"],
        resources:
            n=1,
            time=lambda wildcards, attempt: 12 * 59 * attempt,
            mem_gb_pt=lambda wildcards, attempt: 48 * attempt
        log:
            "results/logs/dbg_filter_variants_f2/{crossing_id}.log"
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


    rule dbg_extract_strict_biallelic_snps_f2:
        input:
            vcf="results/variants/f2/filtered/{crossing_id}.vcf"
        output:
            snps=temp("results/variants/f2/filtered/{crossing_id}.biallelic-snps.raw.vcf"),
            vcf="results/variants/f2/filtered/{crossing_id}.biallelic-snps.vcf"
        params:
            index=config["ref"]["genome"],
            filter=config["debug"]["call_variants_f2"]["biallelic_snps"]["filter"],
            allowed_missing_fraction=config["debug"]["call_variants_f2"]["biallelic_snps"]["allowed_missing_fraction"],
            java_options=config["gatk_options"]["java_options"]
        resources:
            n=1,
            time=lambda wildcards, attempt: 12 * 59 * attempt,
            mem_gb_pt=lambda wildcards, attempt: 48 * attempt
        log:
            "results/logs/dbg_extract_strict_biallelic_snps_f2/{crossing_id}.log"
        conda:
            "../envs/gatk4.yaml"
        shell:
            """
            gatk SelectVariants \
                 -R {params.index} \
                 -V {input.vcf} \
                 --exclude-filtered \
                 --select-type-to-include SNP \
                 --select-type-to-exclude INDEL \
                 --restrict-alleles-to BIALLELIC \
                 --max-nocall-fraction {params.allowed_missing_fraction} \
                 -O {output.snps}  2> {log} ;\
            if [[ "{params.filter}" != "" ]] ; then \
               gatk VariantFiltration \
                    -R {params.index} \
                    -V {output.snps} \
                    --filter-name "biallelic-snps-filter" \
                    --filter-expression "{params.filter}" \
                    -O {output.vcf} 2>> {log} \;
            else \
               cp {output.snps} {output.vcf} 2>> {log} ;\
            fi
            """
