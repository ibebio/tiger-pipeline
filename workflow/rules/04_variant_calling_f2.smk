rule call_variants_f2:
    input:
        bam="results/rmdup/{f2_sample}.rmdup.bam"
        #get_bams_f2 # get all f2 bams for a {crossing_id} wildcard #TODO need to separate each F2 in the list!
    output:
        vcf=temp("results/variants/f2/monomorphic/{crossing_id}.{f2_sample}.monomorphic.vcf"),
    params:
        index=config["ref"]["genome"],
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 24 * attempt
    log:
        "results/logs/call_variants_f2/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        """
        freebayes \
          --fasta-reference {params.index} \
          --report-monomorphic \
          --use-best-n-alleles 1 \
          {input.bam} \
          > {output.vcf} \
        """

rule create_tab_and_gzip_f2:
    input:
        vcf="results/variants/f2/monomorphic/{crossing_id}.{f2_sample}.monomorphic.vcf"
    output:
        vcf=temp("results/variants/f2/monomorphic/{crossing_id}.{f2_sample}.monomorphic.vcf.gz"),
        tab=temp("results/tiger_analysis/F2.{crossing_id}/tab/{f2_sample}.tabbed.txt")
    params:
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 24 * attempt
    log:
        "results/logs/call_variants_f2/{crossing_id}.{f2_sample}.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        """
        bcftools query -f '%CHROM %POS %REF %ALT [ %AO %RO ]\n' {input.vcf} -o {output.tab} \
        ; \
        gzip -c {input.vcf} > {output.vcf}
        """
