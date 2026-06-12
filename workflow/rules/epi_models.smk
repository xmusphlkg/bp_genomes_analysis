# rules/epi_models.smk — T22-T28: Epidemiological enhancement


rule build_ap_exposure_index:
    """T24: Construct aPExposure index (V1 and V2 if PrnInForm available)."""
    input:
        ph_master = config["ph_master"],
        program_metadata = config["ph_program_metadata"],
    output:
        index = f"{WORK}/epi/ap_exposure_index.tsv",
        sensitivity = f"{WORK}/epi/ap_exposure_sensitivity.pdf",
    params:
        version = config["epi"]["exposure_index_version"],
        lambda_range = config["epi"]["exposure_lambda_range"],
        gamma_range = config["epi"]["exposure_gamma_range"],
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/build_ap_exposure_index.py"


rule ipw_prevalence:
    """T25: IPW-corrected prn- prevalence estimates with boundary sensitivity."""
    input:
        manifest = f"{WORK}/manifest/manifest.tsv",
        missingness_model = f"{WORK}/qc/missingness_model.json",
        ph_master = config["ph_master"],
    output:
        prevalence = f"{WORK}/epi/ipw_prevalence.tsv",
        boundary_figure = f"{WORK}/epi/boundary_sensitivity.pdf",
    params:
        weight_truncation = config["epi"]["ipw_weight_truncation"],
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/ipw_prevalence.py"


rule panel_model:
    """T26: Country-level panel association model with FE sensitivity."""
    input:
        exposure = rules.build_ap_exposure_index.output.index,
        prevalence = rules.ipw_prevalence.output.prevalence,
        ph_master = config["ph_master"],
    output:
        results = f"{WORK}/epi/panel_model_results.tsv",
        diagnostics = f"{WORK}/epi/panel_model_diagnostics.pdf",
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/panel_model.py"


rule its_feasibility:
    """T27: Interrupted time series feasibility assessment per country."""
    input:
        prevalence = rules.ipw_prevalence.output.prevalence,
        ph_master = config["ph_master"],
    output:
        report = f"{WORK}/epi/its_feasibility_report.tsv",
    conda:
        "../../config/env/environment_tool.yml"
    script:
        "../lib/its_feasibility.py"
