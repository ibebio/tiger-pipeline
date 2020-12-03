############################################################
# CREATE COMPLETE MARKER FILES
############################################################

localrules: create_complete_markers, create_corrected_markers, create_marker_all

# Extract biallelic SNPs
rule extract_biallelic_snps_parental_complete:
    input:
        vcf="results/variants/parental/filtered_complete/{parental_sample}.vcf"
    output:
        vcf="results/variants/parental/biallelic_snps_complete/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["variant_calling_parental"]["java_options"],
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/extract_biallelic_snps_parental_complete/{parental_sample}.log"
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
             -O {output.vcf}  2> {log}
        """

        
# Extract CHROM and POS into a text file
rule create_complete_markers:
    input:
        vcf="results/variants/parental/biallelic_snps_complete/{parental_sample}.vcf",
    output:
        marker_file="results/markers/complete/{parental_sample}.SNP.biallelic.complete.txt", # REMOVE  SNP.biallelic.complete from name?
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/create_complete_markers/{parental_sample}.log"
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        SnpSift extractFields \
            {input.vcf} \
            CHROM POS > {output.marker_file} 2> {log}
        """

        
############################################################
# CREATE CORRECTED MARKER FILES FROM COMPLETE MARKER FILES
############################################################

# Extract biallelic SNPs
rule extract_biallelic_snps_parental_corrected:
    input:
        vcf="results/variants/parental/filtered_corrected/{parental_sample}.vcf"
    output:
        vcf="results/variants/parental/biallelic_snps_corrected/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["variant_calling_parental"]["java_options"],
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/extract_biallelic_snps/{parental_sample}.log"
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
             -O {output.vcf}  2> {log}
        """

############################################################
# AGGREGATION RULE, collect all marker files/vcfs and switch to a different wildcard
############################################################
# Aggregate all mapped files
rule create_marker_all:
    input:
        parental_biallelic_snps=get_parental_biallelic_snps_corrected_all,
        parental_complete_markers=get_parental_complete_markers_all,
    output:
        flag_vcf=touch("results/variants/parental/biallelic_snps_corrected/biallelic_snps_corrected.done"),
        flag_marker=touch("results/markers/complete/parental_complete_markers.done")



# Intersect
# parent A has Col-0 allelle, parent B has alternative allele
rule intersect_parental_lines:
    input:
        parental_vcfs_done="results/variants/parental/biallelic_snps_corrected/biallelic_snps_corrected.done",
    output:
        bgzip_parent_ref_vcf=temp("results/variants/parental/isec_parent_lines/{crossing_id}.parent_ref.vcf.gz"),
        bgzip_parent_alt_vcf=temp("results/variants/parental/isec_parent_lines/{crossing_id}.parent_alt.vcf.gz"),
        isec_output_dir=temp(directory("results/variants/parental/isec_parent_lines/isec.{crossing_id}")),
        isec_vcf="results/variants/parental/isec_parent_lines/{crossing_id}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["variant_calling_parental"]["java_options"],
        parent_ref_vcf=lambda wildcards: "results/variants/parental/biallelic_snps_corrected/{parent_ref}.vcf".format(parent_ref=[c["parent_ref"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
        parent_alt_vcf=lambda wildcards: "results/variants/parental/biallelic_snps_corrected/{parent_alt}.vcf".format(parent_alt=[c["parent_alt"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/intersect_parental_lines/{crossing_id}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bgzip -c {params.parent_ref_vcf} > {output.bgzip_parent_ref_vcf} 2> {log} ;\
        bcftools index {output.bgzip_parent_ref_vcf} 2> {log} \;
        bgzip -c {params.parent_alt_vcf} > {output.bgzip_parent_alt_vcf} 2> {log} ;\
        bcftools index {output.bgzip_parent_alt_vcf} ;\
        \
        bcftools isec \
            -Ov {output.bgzip_parent_alt_vcf} {output.bgzip_parent_ref_vcf} -p {output.isec_output_dir} ;\
        cp {output.isec_output_dir}/0000.vcf {output.isec_vcf}
        """

rule remove_transposable_elements:
    input:
        vcf="results/variants/parental/isec_parent_lines/{crossing_id}.vcf",
    output:
        vcf="results/variants/parental/isec_parent_lines_noTE/{crossing_id}.vcf",
    params:
        te_bed_file=config["variant_filtering_parental"]["te_bed_file"],
        vcf_prefix="results/variants/parental/isec_parent_lines_noTE/{crossing_id}"
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/remove_transposable_elements/{crossing_id}.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --vcf {input.vcf} --exclude-bed {params.te_bed_file} --recode --out {params.vcf_prefix} ;\
        mv {params.vcf_prefix}.recode.vcf {output.vcf}
        """


# Extract CHROM and POS into a text file
rule create_corrected_markers:
    input:
        vcf="results/variants/parental/isec_parent_lines_noTE/{crossing_id}.vcf"
    output:
        marker_file="results/markers/corrected/{crossing_id}.SNP.private.txt"
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/create_corrected_markers/{crossing_id}.log"
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        SnpSift extractFields \
            {input.vcf} \
            CHROM POS > {output.marker_file} 2> {log}
        """
