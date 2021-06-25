############################################################
# CREATE COMPLETE MARKER FILES
############################################################

localrules: create_complete_markers, create_corrected_markers, create_marker_all, qc_marker_snp_counts

# Extract biallelic SNPs
rule extract_biallelic_snps_parental_complete:
    input:
        vcf="results/variants/parental/filtered_complete/{parental_sample}.vcf"
    output:
        vcf="results/variants/parental/biallelic_snps_complete/{parental_sample}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["gatk_options"]["java_options"],
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
        java_options=config["gatk_options"]["java_options"],
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
# The ref and alt parents will either be determined automatically, based on which of the parents, a or b, has the
# higher number of SNPs, or manually. In the manual case, parent_a is set to ref and parent_b to alt.
# parent A has Col-0 allele, parent B has alternative allele
rule intersect_parental_lines:
    input:
        parental_vcfs_done="results/variants/parental/biallelic_snps_corrected/biallelic_snps_corrected.done",
    output:
        bgzip_parent_ref_vcf=temp("results/variants/parental/isec_parent_lines/{crossing_id}.parent_ref.vcf.gz"),
        bgzip_parent_src_vcf=temp("results/variants/parental/isec_parent_lines/{crossing_id}.parent_src.vcf.gz"),
        isec_output_dir=temp(directory("results/variants/parental/isec_parent_lines/isec.{crossing_id}")),
        isec_output_file=temp("results/variants/parental/isec_parent_lines/isec.{crossing_id}/0000.vcf"),
        isec_vcf="results/variants/parental/isec_parent_lines/{crossing_id}.vcf",
    params:
        index=config["ref"]["genome"],
        java_options=config["gatk_options"]["java_options"],
        parent_a_vcf=lambda wildcards: "results/variants/parental/biallelic_snps_corrected/{parent_a}.vcf".format(parent_a=[c["parent_a"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
        parent_b_vcf=lambda wildcards: "results/variants/parental/biallelic_snps_corrected/{parent_b}.vcf".format(parent_b=[c["parent_b"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
        auto_order=lambda wildcards: "{auto_order}".format(auto_order=[c["auto_order"] for c in config["crossings"] if c["id"] == wildcards.crossing_id][0]),
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
        # Select marker ref and marker source
        if [[ "{params.auto_order}" == "yes" ]] ; then
           echo "Try to automatically determine marker source and marker ref parental lines" 2> {log} ;\
           a_snp_count=$(cat {params.parent_a_vcf} |grep -v '^#' |wc -l) ;\
           b_snp_count=$(cat {params.parent_b_vcf} |grep -v '^#' |wc -l) ;\
           if [[ $a_snp_count -gt $b_snp_count ]] ; then \
             parent_src_vcf={params.parent_a_vcf} ;\
             parent_ref_vcf={params.parent_b_vcf} ;\
             echo "Selected parent a as marker source with $a_snp_count SNPs" 2>> {log} ;\
             echo "Selected parent b as marker ref with $b_snp_count SNPs" 2>> {log} ;\
           else \
             parent_src_vcf={params.parent_b_vcf} ;\
             parent_ref_vcf={params.parent_a_vcf} ;\
             echo "Selected parent b as marker source with $b_snp_count SNPs" 2>> {log} ;\
             echo "Selected parent a as marker ref with $a_snp_count SNPs" 2>> {log} ;\
           fi ;\
        else \
           echo "Using manual ref/src order: parent a is ref, parent b is src."
           parent_src_vcf={params.parent_b_vcf} ;\
           parent_ref_vcf={params.parent_a_vcf} ;\
        fi ;\
        \
        # Prepare parent a and b input files
        bgzip -c $parent_ref_vcf > {output.bgzip_parent_ref_vcf} 2>> {log} ;\
        bcftools index {output.bgzip_parent_ref_vcf} 2>> {log} \;
        bgzip -c $parent_src_vcf > {output.bgzip_parent_src_vcf} 2>> {log} ;\
        bcftools index {output.bgzip_parent_src_vcf} ;\
        \
        bcftools isec \
            -Ov {output.bgzip_parent_src_vcf} {output.bgzip_parent_ref_vcf} -p {output.isec_output_dir} ;\
        cp {output.isec_output_dir}/0000.vcf {output.isec_vcf}
        """


# Create figure with the number of SNPs
rule qc_marker_snp_counts:
    input:
        isec_output="results/variants/parental/isec_parent_lines/isec.{crossing_id}/0000.vcf",
        isec_output_dir=directory("results/variants/parental/isec_parent_lines/isec.{crossing_id}"),
        # corrected_marker_vcf="results/variants/parental/isec_parent_lines_noTE/{crossing_id}.vcf"
    output:
        parental_snps_fig="results/qc/biallelic_parental_snp_counts_{crossing_id}_mqc.png",
    params:
        crossing_id="{crossing_id}"
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 1 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 24 * attempt
    log:
        "results/logs/qc_marker_snp_counts/{crossing_id}.log"
    conda:
        "../envs/qc.yaml"
    script:
        # rename to plot_marker_snp_counts
        "../scripts/plot_parental_biallelic_snp_counts.py"

        
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


rule remove_telomeric_regions:
    input:
        vcf="results/variants/parental/isec_parent_lines_noTE/{crossing_id}.vcf",
    output:
        vcf="results/variants/parental/isec_parent_lines_noTEnoTLR/{crossing_id}.vcf",
    params:
        tlr_bed_file=config["variant_filtering_parental"]["tlr_bed_file"],
        vcf_prefix="results/variants/parental/isec_parent_lines_noTEnoTLR/{crossing_id}"
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt
    log:
        "results/logs/remove_telomeric_regions/{crossing_id}.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --vcf {input.vcf} --exclude-bed {params.tlr_bed_file} --recode --out {params.vcf_prefix} ;\
        mv {params.vcf_prefix}.recode.vcf {output.vcf}
        """


# Extract CHROM and POS into a text file
rule create_corrected_markers:
    input:
        vcf="results/variants/parental/isec_parent_lines_noTEnoTLR/{crossing_id}.vcf"
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
