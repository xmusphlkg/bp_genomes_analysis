# rules/snippy.smk — T11: Reference-based SNP calling with Snippy


rule snippy_per_sample_reads:
    """Run Snippy per sample (reads mode — preferred)."""
    input:
        r1 = f"{WORK}/reads_clean/{{sample}}_1.fastq.gz",
        r2 = f"{WORK}/reads_clean/{{sample}}_2.fastq.gz",
        ref = REF_GENOME,
    output:
        outdir = directory(f"{WORK}/snippy/{{sample}}"),
        vcf = f"{WORK}/snippy/{{sample}}/snps.vcf",
    params:
        mincov = config["snippy"]["mincov"],
        minfrac = config["snippy"]["minfrac"],
        mapqual = config["snippy"]["mapqual"],
        basequal = config["snippy"]["basequal"],
    threads: config["snippy"]["cpus"]
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        snippy \
            --cpus {threads} \
            --outdir {output.outdir} \
            --ref {input.ref} \
            --R1 {input.r1} --R2 {input.r2} \
            --mincov {params.mincov} \
            --minfrac {params.minfrac} \
            --mapqual {params.mapqual} \
            --basequal {params.basequal} \
            --force
        """


rule snippy_per_sample_contigs:
    """Run Snippy per sample (contig mode — fallback for samples without reads)."""
    input:
        contigs = f"{config['genome_fasta_dir']}/{{sample}}.fna",
        ref = REF_GENOME,
    output:
        outdir = directory(f"{WORK}/snippy_ctg/{{sample}}"),
        vcf = f"{WORK}/snippy_ctg/{{sample}}/snps.vcf",
    threads: config["snippy"]["cpus"]
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        snippy \
            --cpus {threads} \
            --outdir {output.outdir} \
            --ref {input.ref} \
            --ctgs {input.contigs} \
            --force
        """


rule snippy_core:
    """Build core SNP alignment from all per-sample Snippy outputs."""
    input:
        ref = REF_GENOME,
        snippy_dirs = expand(f"{WORK}/snippy/{{s}}", s=SAMPLES_WITH_READS)
                    + expand(f"{WORK}/snippy_ctg/{{s}}", s=SAMPLES_CONTIGS_ONLY),
    output:
        full_aln = f"{WORK}/phylo/core.full.aln",
        aln = f"{WORK}/phylo/core.aln",
        txt = f"{WORK}/phylo/core.txt",
    conda:
        "../../config/env/environment_tool.yml"
    shell:
        """
        snippy-core \
            --ref {input.ref} \
            --prefix {WORK}/phylo/core \
            {input.snippy_dirs}
        """
