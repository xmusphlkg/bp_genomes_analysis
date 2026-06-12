# rules/asr.smk — T15/T16/T17: Ancestral state reconstruction + origin events


rule asr_parsimony:
    """Fitch parsimony ancestral state reconstruction."""
    input:
        tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        manifest = f"{WORK}/manifest/manifest.tsv",
    output:
        states = f"{WORK}/asr/parsimony_states.tsv",
        transitions = f"{WORK}/asr/parsimony_transitions.tsv",
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/asr_parsimony.py"


rule asr_pastml:
    """PastML likelihood ancestral state reconstruction."""
    input:
        tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        manifest = f"{WORK}/manifest/manifest.tsv",
    output:
        states = f"{WORK}/asr/pastml_states.tsv",
        html = f"{WORK}/asr/pastml_visualization.html",
    params:
        model = config["asr"]["pastml_model"],
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/asr_pastml.py"


rule asr_multi_scenario:
    """Run ASR under multiple missing-data scenarios (mask / conservative / impute)."""
    input:
        tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        manifest = f"{WORK}/manifest/manifest.tsv",
        missingness_model = f"{WORK}/qc/missingness_model.json",
    output:
        summary = f"{WORK}/asr/multi_scenario_summary.tsv",
    params:
        scenarios = config["asr"]["missing_scenarios"],
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/asr_multi_scenario.py"


rule origin_events:
    """T16: Identify independent origin events with evidence packages."""
    input:
        pars_states = rules.asr_parsimony.output.states,
        pars_transitions = rules.asr_parsimony.output.transitions,
        pastml_states = rules.asr_pastml.output.states,
        tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        manifest = f"{WORK}/manifest/manifest.tsv",
        mechanism_calls = config["step4_mechanism_calls"],
    output:
        events = f"{WORK}/asr/origin_events.tsv",
        event_trees = directory(f"{WORK}/asr/event_subtrees/"),
    params:
        min_bootstrap = config["asr"]["min_bootstrap_threshold"],
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/origin_events.py"


rule robustness_matrix:
    """T17: Multi-method robustness testing."""
    input:
        iq_tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        rax_tree = f"{WORK}/phylo/raxmlng/ml_tree.raxml.bestTree",
        multi_scenario = rules.asr_multi_scenario.output.summary,
        origin_events = rules.origin_events.output.events,
    output:
        matrix = f"{WORK}/asr/robustness_matrix.tsv",
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/robustness_matrix.py"
