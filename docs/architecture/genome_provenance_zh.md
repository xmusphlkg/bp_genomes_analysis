# 基因组来源与路径说明

## 为什么要区分 `public genome paths` 和 `raw-read genome paths`

它们都可能追溯到 NCBI / ENA / SRA，但不是同一种数据产物。

- `public genome paths`
  - 指向的是公共数据库已经发布好的组装 FASTA
  - 典型 accession 是 `GCA_*` / `GCF_*`
  - 当前 lookup 表是 `modules/step2_typing/outputs/bp_genome_paths_qc.tsv`
  - 这张表由 `modules/step1_ingest/outputs/bp_public_genome_qc_manifest.tsv` 派生，并把路径固定到 `pertussis_data/bp_genomes_qc/assemblies/`
  - `pertussis_data/bp_genomes_qc/assemblies/` 是 public assembly 的正式本地落点，不是临时兜底目录
  - 这些 FASTA 是“直接下载公共组装”，再由正式 registry 统一指向

- `raw-read genome paths`
  - 指向的是我们从公共 raw reads 下载 FASTQ 后，在本地重新组装得到的 FASTA
  - 典型 accession 是仓库内部生成的 `RRASM_*`
  - 当前主表是 `modules/step1_ingest/outputs/bp_raw_read_step3_genome_paths.tsv`
  - 这些 FASTA 是“公共 reads -> 本地 assembly”

## 为什么不能混成同一张 path 表

必须分开，原因有 3 个：

1. provenance 不同  
   一个是公共数据库直接发布的 assembly，一个是我们本地重建的 assembly。

2. 去重逻辑不同  
   同一个样本可能同时有公共 assembly 和 raw reads；如果不分开，后面很容易把“同一样本的两个不同 genome product”混成一条。

3. 分析解释不同  
   `public_read_rescue` 表示我们保留的是公共 assembly，只是 reads linkage 被补全了；  
   `public_raw_read_assembly` 才表示最终进入分析的是本地 raw-read assembly。

## 现在新增的统一总表

为了避免来回查多张表，现在新增：

- `outputs/workflow/manifest/genome_catalog.tsv`

这张表以 canonical `manifest.tsv` 为底表，额外增加：

- `genome_path_registry`
- `genome_file_class`
- `primary_fasta_accession`
- `primary_fasta_path`
- `primary_fasta_status`
- `primary_fasta_note`

其中：

- `genome_file_class = direct_public_assembly_fasta`
  - 说明最终 FASTA 来自公共 assembly

- `genome_file_class = local_raw_read_assembly_fasta`
  - 说明最终 FASTA 来自本地 raw-read assembly

- `genome_path_registry = unresolved`
  - 说明 manifest 里有这个 genome，但现有 path registry 还没有完整回填到这一条
  - 这类记录值得继续补 provenance，但不会改变它在 manifest 里的样本级元数据

- `genome_path_registry = bp_step2_public_genome_paths_qc`
  - 说明 public assembly 的实际 FASTA 来自正式的 public path registry
  - 这张 registry 的路径最终仍指向 `pertussis_data/bp_genomes_qc/assemblies/`
  - `primary_fasta_note` 会保留 `public_registry_lookup=hit|miss`，用于记录旧的 lookup 表是否命中

- `genome_path_registry = bp_genomes_qc_assemblies`
  - 这是兼容性回退分支，仅在正式 registry 仍缺条目时才会出现
  - 当前正式 registry 已覆盖全部 public 记录，正常重建时不应再看到它

- `primary_fasta_note`
  - 对 public assembly 行，通常会带上 `resolved_from=pertussis_data/bp_genomes_qc/assemblies`
  - 对命中旧 lookup 表但实际路径仍需回填的行，还会显示 `public_registry_lookup=hit` 和 `public_registry_status=...`

## 与 manuscript 用表的区别

- `manuscript/figure_data/project_genome_metadata_manifest.tsv`
  - 是 submission-facing 的简化样本表
  - 适合论文数据说明

- `outputs/workflow/manifest/genome_catalog.tsv`
  - 是技术审计用全量目录
  - 适合排查 provenance、路径、去重和分析输入来源
