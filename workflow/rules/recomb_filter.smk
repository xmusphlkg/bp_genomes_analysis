# rules/recomb_filter.smk — T12: Recombination filtering (Gubbins + ClonalFrameML)


rule gubbins:
    """Identify recombination and build non-recombinant ML tree."""
    input:
        aln = f"{WORK}/phylo/core.full.aln",
    output:
        filtered_aln = f"{WORK}/phylo/gubbins/core.filtered_polymorphic_sites.fasta",
        tree = f"{WORK}/phylo/gubbins/core.final_tree.tre",
        recomb_gff = f"{WORK}/phylo/gubbins/core.recombination_predictions.gff",
        summary = f"{WORK}/phylo/gubbins/core.per_branch_statistics.csv",
    params:
        iterations = config["recomb_filter"]["gubbins_iterations"],
        prefix = f"{WORK}/phylo/gubbins/core",
    threads: 8
    conda:
        "../../config/env/environment_python.yml"
    shell:
        """
        run_gubbins.py \
            --prefix {params.prefix} \
            --threads {threads} \
            --iterations {params.iterations} \
            --tree-builder iqtree \
            --mar \
            --seq-recon iqtree \
            {input.aln}
        """


rule clonalframeml:
    """Cross-validate recombination inference with ClonalFrameML."""
    input:
        tree = rules.gubbins.output.tree,
        aln = f"{WORK}/phylo/core.full.aln",
    output:
        cfml_tree = f"{WORK}/phylo/clonalframeml/core.labelled_tree.newick",
        em_results = f"{WORK}/phylo/clonalframeml/core.em.txt",
        importation = f"{WORK}/phylo/clonalframeml/core.importation_status.txt",
    params:
        prefix = f"{WORK}/phylo/clonalframeml/core",
        emsim = config["recomb_filter"]["clonalframeml_emsim"],
    conda:
        "../../config/env/environment_python.yml"
    shell:
        """
        ClonalFrameML \
            {input.tree} \
            {input.aln} \
            {params.prefix} \
            -emsim {params.emsim}
        """


rule mask_recombination:
    """Generate recomb-filtered alignment for downstream ML tree building."""
    input:
        full_aln = f"{WORK}/phylo/core.full.aln",
        gubbins_gff = rules.gubbins.output.recomb_gff,
    output:
        filtered_aln = f"{WORK}/phylo/recomb_filtered.aln",
    conda:
        "../../config/env/environment_python.yml"
    script:
        "../lib/mask_recombination.py"
