# rules/ml_tree.smk — T13: ML phylogeny (IQ-TREE2 + RAxML-NG cross-validation)


rule iqtree2:
    """Maximum likelihood tree with IQ-TREE2 (primary)."""
    input:
        aln = f"{WORK}/phylo/recomb_filtered.aln",
    output:
        tree = f"{WORK}/phylo/iqtree2/ml_tree.treefile",
        log = f"{WORK}/phylo/iqtree2/ml_tree.log",
        iqtree = f"{WORK}/phylo/iqtree2/ml_tree.iqtree",
    params:
        model = config["ml_tree"]["iqtree2_model"],
        bootstrap = config["ml_tree"]["iqtree2_bootstrap"],
        prefix = f"{WORK}/phylo/iqtree2/ml_tree",
    threads: 16
    conda:
        "../../config/env/environment_python.yml"
    shell:
        """
        iqtree2 \
            -s {input.aln} \
            -m {params.model} \
            -bb {params.bootstrap} \
            -nt {threads} \
            --prefix {params.prefix} \
            -redo
        """


rule raxmlng:
    """Maximum likelihood tree with RAxML-NG (cross-validation)."""
    input:
        aln = f"{WORK}/phylo/recomb_filtered.aln",
    output:
        tree = f"{WORK}/phylo/raxmlng/ml_tree.raxml.bestTree",
        support = f"{WORK}/phylo/raxmlng/ml_tree.raxml.support",
        log = f"{WORK}/phylo/raxmlng/ml_tree.raxml.log",
    params:
        model = config["ml_tree"]["raxmlng_model"],
        bootstrap = config["ml_tree"]["raxmlng_bootstrap"],
        prefix = f"{WORK}/phylo/raxmlng/ml_tree",
    threads: 16
    conda:
        "../../config/env/environment_python.yml"
    shell:
        """
        raxml-ng --all \
            --msa {input.aln} \
            --model {params.model} \
            --bs-trees {params.bootstrap} \
            --threads {threads} \
            --prefix {params.prefix} \
            --redo
        """


rule compare_trees:
    """Compare IQ-TREE2 and RAxML-NG topologies for consistency."""
    input:
        iq_tree = rules.iqtree2.output.tree,
        rax_tree = rules.raxmlng.output.tree,
    output:
        comparison = f"{WORK}/phylo/tree_comparison_report.json",
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/compare_trees.py"
