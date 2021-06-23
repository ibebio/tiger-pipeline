rule collect_allelic_counts_f2:
    input:
        bam="results/rmdup/{f2_sample}.rmdup.bam",
        marker_complete_ref_vcf=lambda wildcards: "results/variants/parental/biallelic_snps_complete/{parent_a}.vcf".format(parent_a=[c["parent_a"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0])
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
            -L {input.marker_complete_ref_vcf} \
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




# rule create_tab_and_gzip_f2:
#     input:
#         vcf="results/variants/f2/allelic_counts/{crossing_id}.{f2_sample}."
#     output:
#         vcf=temp("results/variants/f2/monomorphic/{crossing_id}.{f2_sample}.monomorphic.vcf.gz"),
#         tab=temp("results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tabbed.txt")
#     params:
#     resources:
#         n=1,
#         time=lambda wildcards, attempt: 12 * 59 * attempt,
#         mem_gb_pt=lambda wildcards, attempt: 24 * attempt
#     log:
#         "results/logs/call_variants_f2/{crossing_id}.{f2_sample}.log"
#     conda:
#         "../envs/freebayes.yaml"
#     shell:
#         """
#         bcftools query -f '%CHROM %POS %REF %ALT [ %AO %RO ]\n' {input.vcf} -o {output.tab} \
#         ; \
#         gzip -c {input.vcf} > {output.vcf}
#         """
