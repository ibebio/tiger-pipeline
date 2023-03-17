############################################################
# 1st STEP: CREATE .BAM FILES W/O DUPLICATES
# FOR _ALL_ SAMPLES INVOLVED, PARENTAL AND F2
############################################################

localrules: map_all,remove_duplicates

# Trim reads (PE-only)
rule trim_reads:
    priority: 10
    input:
        unpack(get_fastq) # TODO ALL FASTQ'S, parental and f2
    output:
        fq1=temp("results/trimmed/{sample}.R1.trimmed.fastq.gz"),
        fq2=temp("results/trimmed/{sample}.R2.trimmed.fastq.gz"),
    params:
        output_dir="results/trimmed",
        basename="{sample}"
    threads: 4
    resources:
        n=4,
        time=lambda wildcards, attempt: 2 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 6 * attempt,
    log:
        "results/logs/trim_reads/{sample}.log"
    conda:
        "../envs/trim_reads.yaml"
    shell:
        """
        trim_galore \
        --cores {threads} \
        --gzip \
        --paired \
        --output_dir {params.output_dir} \
        --basename {params.basename} \
        {input} 2> {log} ; \
        mv {params.output_dir}/{params.basename}_val_1.fq.gz {output.fq1} ; \
        mv {params.output_dir}/{params.basename}_val_2.fq.gz {output.fq2}
        """

#(bio36) ibezrukov2@taco:/ebio/abt6_projects9/At_dFLC/data/dFLC_F2mapping_AS$ fastqc -o ~/fastq_test/ p9925-WT-001.merged2.R1.fastq.gz
#p9925-WT-001.merged2.R1_fastqc.html  p9925-WT-001.merged2.R1_fastqc.zip
        
        
rule qc_trim_reads:
    input:
        fq1="results/trimmed/{sample}.R1.trimmed.fastq.gz",
        fq2="results/trimmed/{sample}.R2.trimmed.fastq.gz",
    output:
        qc_zip1="results/qc/trimmed_fastqc/{sample}.R1.trimmed_fastqc.zip",
        qc_zip2="results/qc/trimmed_fastqc/{sample}.R2.trimmed_fastqc.zip",
        qc_html1="results/qc/trimmed_fastqc/{sample}.R1.trimmed_fastqc.html",
        qc_html2="results/qc/trimmed_fastqc/{sample}.R2.trimmed_fastqc.html",        
    threads: 2
    params:
        outdir="results/qc/trimmed_fastqc/"
    resources:
        n=2,
        time=lambda wildcards, attempt: 1 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 16 * attempt,
    log:
        "results/logs/qc_trim_reads/{sample}.log"
    conda:
        "../envs/qc.yaml"
                        
    shell:
         "fastqc -t {threads} -o {params.outdir} {input.fq1} {input.fq2} 2> {log}"


# Align to reference
rule map_to_reference:
    priority: 9
    input:
        fq1="results/trimmed/{sample}.R1.trimmed.fastq.gz",
        fq2="results/trimmed/{sample}.R2.trimmed.fastq.gz",
        qc_zip1="results/qc/trimmed_fastqc/{sample}.R1.trimmed_fastqc.zip", # The QC files are not used, but should be
        qc_zip2="results/qc/trimmed_fastqc/{sample}.R2.trimmed_fastqc.zip", # created before running this rule.
    output:
        temp("results/mapped/{sample}.sorted.bam")
    params:
        index=config["ref"]["genome"],
        platform=config["read_group"]["platform"],
        library=config["read_group"]["library"]
        # read_group: geat_read_group,
    threads: 10
    resources:
        n=10,
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 4 * attempt,
    log:
        "results/logs/map_to_reference/{sample}.log"
    conda:
        "../envs/global.yaml"
                        
    shell:
        "workflow/scripts/map_to_reference.sh {threads} {params.platform} {wildcards.sample} {params.index} {output} {input} {params.library} 2> {log}"


rule remove_duplicates:
    priority: 8
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam="results/rmdup/{sample}.rmdup.bam",
        metrics="results/qc/rmdup/{sample}.metrics.txt"
    threads: 2
    resources:
        time=lambda wildcards, attempt: 3 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 32 * attempt,
    log:
        "results/logs/remove_duplicates/{sample}.log"
    conda:
        "../envs/rmdup.yaml"
    shell:
        """
        picard MarkDuplicates \
           REMOVE_DUPLICATES=true \
           VALIDATION_STRINGENCY=LENIENT \
           INPUT={input} \
           METRICS_FILE={output.metrics} \
           OUTPUT={output.bam} \
        ; \
        samtools index \
           {output.bam} 2> {log}
        """

rule qc_bam_qualimap:
    input:
        bam="results/rmdup/{sample}.rmdup.bam",
    output:
        qc_dir=directory("results/qc/qualimap/{sample}"),
        qc_genome_results="results/qc/qualimap/{sample}/genome_results.txt",
    resources:
        time=lambda wildcards, attempt: 1 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 38 * attempt,
    log:
        "results/logs/qc_bam_qualimap/{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        qualimap bamqc --java-mem-size=16G \
                       -bam {input.bam} \
                       -outdir {output.qc_dir} \
        2> {log}
        """

# Aggregate all mapped files
rule map_all:
    input:
        bams=get_all_mapped_files,
        qcs=expand("results/qc/qualimap/{sample}/genome_results.txt", sample=samples.index)
        
    output:
        flag=touch("results/rmdup/rmdup.done")
