# rules/reads_qc.smk — T07: Reads QC + preprocessing


rule fastp_qc:
    """Per-sample read QC and adapter trimming with fastp."""
    input:
        r1 = f"{WORK}/reads/{{sample}}_1.fastq.gz",
        r2 = f"{WORK}/reads/{{sample}}_2.fastq.gz",
    output:
        r1 = f"{WORK}/reads_clean/{{sample}}_1.fastq.gz",
        r2 = f"{WORK}/reads_clean/{{sample}}_2.fastq.gz",
        json = f"{WORK}/qc/fastp/{{sample}}.fastp.json",
        html = f"{WORK}/qc/fastp/{{sample}}.fastp.html",
    threads: 2
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --json {output.json} --html {output.html} \
            --thread {threads} \
            --qualified_quality_phred 15 \
            --length_required 50 \
            --detect_adapter_for_pe
        """


rule kraken2_screen:
    """Contamination screening with Kraken2."""
    input:
        r1 = rules.fastp_qc.output.r1,
        r2 = rules.fastp_qc.output.r2,
    output:
        report = f"{WORK}/qc/kraken2/{{sample}}.kraken2.report",
        output_seqs = f"{WORK}/qc/kraken2/{{sample}}.kraken2.out",
    params:
        db = config.get("kraken2_db", "references/kraken2_db"),
    threads: 4
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --paired {input.r1} {input.r2} \
            --report {output.report} \
            --output {output.output_seqs} \
            --threads {threads} \
            --confidence 0.2
        """


rule multiqc_aggregate:
    """Aggregate all per-sample QC into a single MultiQC report."""
    input:
        fastp = expand(f"{WORK}/qc/fastp/{{s}}.fastp.json", s=SAMPLES_WITH_READS),
        kraken = expand(f"{WORK}/qc/kraken2/{{s}}.kraken2.report", s=SAMPLES_WITH_READS),
    output:
        html = f"{WORK}/qc/multiqc_report.html",
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        multiqc {WORK}/qc/ --outdir {WORK}/qc/ --filename multiqc_report --force
        """
