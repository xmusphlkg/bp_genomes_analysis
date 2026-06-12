# rules/assembly_qc.smk — T08: Assembly QC (all 2,247 standardized)


rule quast_qc:
    """Assembly QC with QUAST against Tohama I reference."""
    input:
        assembly = f"{config['genome_fasta_dir']}/{{sample}}.fna",
        reference = REF_GENOME,
    output:
        report = f"{WORK}/qc/quast/{{sample}}/report.tsv",
    threads: 2
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        quast.py {input.assembly} \
            -r {input.reference} \
            -o {WORK}/qc/quast/{wildcards.sample} \
            --threads {threads} \
            --min-contig 500
        """


rule aggregate_assembly_qc:
    """Merge per-sample QUAST + contamination results into one table."""
    input:
        quast = expand(f"{WORK}/qc/quast/{{s}}/report.tsv", s=ALL_SAMPLES),
    output:
        table = f"{WORK}/qc/assembly_qc_report.tsv",
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/aggregate_assembly_qc.py"


rule missingness_model:
    """T09/T25: Model determinants of prn interpretability (selection bias diagnostics)."""
    input:
        manifest = f"{WORK}/manifest/manifest.tsv",
        assembly_qc = rules.aggregate_assembly_qc.output.table,
    output:
        model = f"{WORK}/qc/missingness_model.json",
        report = f"{WORK}/qc/missingness_diagnostics.html",
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/missingness_model.py"
