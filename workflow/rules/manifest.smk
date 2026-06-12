# rules/manifest.smk — T02: Unified manifest (single source of truth)


rule build_manifest:
    """Unify sample metadata from Step1/Step2/Step4/SupTable1 into a single manifest."""
    input:
        step1_meta   = config["step1_metadata"],
        step2_qc     = config["step2_qc_table"],
        step2_typing = config["step2_typing_manifest"],
        step4_mech   = config["step4_mechanism_calls"],
        supp_table1  = config["supp_table1"],
    output:
        manifest     = f"{WORK}/manifest/manifest.tsv",
        report       = f"{WORK}/manifest/manifest_build_report.json",
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/build_analysis_manifest.py"


rule trace_reads_availability:
    """T04: Reads availability for each sample via BioSample/BioProject."""
    input:
        manifest = rules.build_manifest.output.manifest,
    output:
        runs_table = f"{WORK}/manifest/runs.tsv",
        reads_report = f"{WORK}/checkpoints/reads_availability_report.json",
    params:
        min_reads_pct = config["reads_availability_min_pct"],
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../../modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py"
