# rules/is_detection.smk — thin delegate into Step4-local PRN read validation assets


rule step4_prn_read_validation:
    """Delegate Step4 PRN read-based IS validation to stage-local scripts and outputs."""
    input:
        subset = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_validation_subset.tsv"),
        mechanism_calls = config["step4_mechanism_calls"],
        is_reference = config.get(
            "is_reference",
            str(PROJECT_DATA_ROOT / "step4_prn_validation" / "references" / "is_elements" / "bp_is_reference.fasta"),
        ),
    output:
        validation = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation.tsv"),
        evidence = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation_is_calls.tsv"),
        tsd = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_read_validation_tsd.tsv"),
        summary = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_validation_summary.tsv"),
        hotspot = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_is_hotspot_results.tsv"),
        hotspot_figure = str(PROJECT_DATA_ROOT / "step4_prn_validation" / "outputs" / "bp_prn_is_hotspot_density.pdf"),
    params:
        batch_label = "snakemake",
        step4_root = str(PROJECT_DATA_ROOT / "step4_prn_validation"),
        hotspot_permutations = config["is_detection"]["hotspot_permutations"],
        panisa_min_support = config["is_detection"]["panisa_min_support"],
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        bash ../../modules/step4_prn_validation/bin/step4_03e_run_is_read_validation.sh \
            --batch-label {params.batch_label} \
            --min-support {params.panisa_min_support}
        python ../../modules/step4_prn_validation/bin/step4_03_validate_prn_with_reads.py \
            --batch {params.step4_root}/work/read_validation/{params.batch_label}/bp_prn_read_validation_batch.tsv \
            --batch-label {params.batch_label} \
            --is-work-root {params.step4_root}/work/read_validation/{params.batch_label} \
            --out {output.validation} \
            --evidence-out {output.evidence} \
            --tsd-out {output.tsd}
        python ../../modules/step4_prn_validation/bin/step4_04_summarize_prn_validation.py \
            --validation {output.validation} \
            --mechanism-calls {input.mechanism_calls} \
            --out {output.summary}
        python ../../modules/step4_prn_validation/bin/step4_03f_hotspot_test.py \
            --evidence {output.evidence} \
            --n-permutations {params.hotspot_permutations} \
            --out {output.hotspot} \
            --plot {output.hotspot_figure}
        """
