# Step6 references (marker FASTAs)

This folder is intentionally **empty** in git: you should add your own marker query sequences.

Step6 scripts expect nucleotide FASTA files for the markers you want to scan.

## Where to put files

- Put marker query FASTAs under `references/markers/`.

Example filenames (recommended):

- `references/markers/prn_maker.fasta`
- `references/markers/ptxP_promoter.fasta`
- `references/markers/fim2.fasta`
- `references/markers/fim3.fasta`

For macrolide resistance screening (23S rRNA A2047G), also provide:

- `references/23S_rRNA.fasta`

Notes:

- These queries should be **B. pertussis** sequences in the correct gene/promoter orientation.
- For allele tracking, Step6 will extract the best-hit subject sequence and assign an allele ID by `md5` hash.
- For 23S, Step6 maps the reference position (default 2047, 1-based) to the aligned subject base.
