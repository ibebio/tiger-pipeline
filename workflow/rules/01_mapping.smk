############################################################
# 1st STEP: CREATE .BAM FILES W/O DUPLICATES
# FOR _ALL_ SAMPLES INVOLVED, PARENTAL AND F2
############################################################

localrules: map_all

# Trim reads (PE-only)
rule trim_reads:
    input:
        unpack(get_fastq) # TODO ALL FASTQ'S, parental and f2
    output:
        fq1=temp("results/trimmed/{sample}.R1.trimmed.fastq.gz"),
        fq2=temp("results/trimmed/{sample}.R2.trimmed.fastq.gz"),
    params:
        adapter=config["read_trimming"]["adapter"]
    threads: 1
    resources:
        n=1,
        time=lambda wildcards, attempt: 2 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt,
    log:
        "results/logs/trim_reads/{sample}.log"
    conda:
        "../envs/trim_reads.yaml"
    shell:
        """
        cutadapt \
        -a {params.adapter} \
        -o {output.fq1} \
        -p {output.fq2} \
        {input}
        2> {log}
        """


# Align to reference
rule map_to_reference:
    input:
        fq1="results/trimmed/{sample}.R1.trimmed.fastq.gz",
        fq2="results/trimmed/{sample}.R2.trimmed.fastq.gz",
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
        mem_gb_pt=lambda wildcards, attempt: 2 * attempt,
    log:
        "results/logs/map_to_reference/{sample}.log"
    conda:
        "../envs/global.yaml"
                        
    shell:
        "workflow/scripts/map_to_reference.sh {threads} {params.platform} {wildcards.sample} {params.index} {output} {input} {params.library} 2> {log}"

        
rule remove_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam="results/rmdup/{sample}.rmdup.bam",
        metrics="results/rmdup/{sample}/metrics.txt"
    resources:
        time=lambda wildcards, attempt: 12 * 59 * attempt,
        mem_gb_pt=lambda wildcards, attempt: 12 * attempt,
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

# Aggregate all mapped files
rule map_all:
    input:
        get_all_mapped_files
    output:
        flag=touch("results/rmdup/rmdup.done")
