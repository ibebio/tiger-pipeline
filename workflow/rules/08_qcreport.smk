localrules: multiqc_report

rule multiqc_report:
    input:
        qc_fastqc_zip1=expand("results/qc/trimmed_fastqc/{sample}.R1.trimmed_fastqc.zip", sample=samples.index),
        qc_fastqc_zip2=expand("results/qc/trimmed_fastqc/{sample}.R2.trimmed_fastqc.zip", sample=samples.index),
        qc_bam_qualimap=expand("results/qc/qualimap/{sample}/genome_results.txt", sample=samples.index),
        qc_marker_snp_counts=expand("results/qc/biallelic_parental_snp_counts_{crossing_id}_mqc.png", crossing_id=[c["id"] for c in config["crossings"]])
    output:
        multiqc_report="results/report/qc.html"
    params:
        qc_root_dir="results/qc/"
    log:
        "results/logs/multiqc_report/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.qc_root_dir} -n {output.multiqc_report} 2> {log}
        """

    
