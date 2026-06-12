# bp_step1 script layout

`step1_ingest/bin` 已按职责拆成 3 个子目录：

- `core/`
  `01_fetch_ncbi_report.py` 到 `05_verify_genome_extract.py`
  负责 NCBI Assembly 元数据抓取、清洗、下载与基础校验。
- `manifest/`
  `06_build_public_manifest.py` 到 `09_build_analysis_cohorts.py`
  负责 public genome manifest、raw-read 链接恢复、去重和 cohort 构建。
- `raw_reads/`
  `10_build_download_plan.py` 到 `19_merge_qc_passed_genomes.py`，以及 `20_build_genome_database.py` 到 `27_run_snippy_batch.sh`
  负责 SRA/ENA raw-read 扩增、分片、分布式下载装配、环境预检、external gap-fill、装配后 QC、合并、统一 genome database 的补全与完整性维护，以及 manifest 驱动的 reads tracing / reads planning / assembly QC / Snippy 批处理。

常用入口：

- `bash workflow/bin/run_full_workflow.sh --dry-run`
- `python3 modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py`
- `python3 modules/step1_ingest/bin/raw_reads/10_build_download_plan.py`
- `python3 modules/step1_ingest/bin/raw_reads/11_split_download_plan.py`
- `python3 modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py`
- `python3 modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py`
- `python3 modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py`
- `python3 modules/step1_ingest/bin/raw_reads/20_build_genome_database.py`
- `bash modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh`
- `python3 modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py`
- `bash modules/step1_ingest/bin/raw_reads/23_retry_missing_assemblies.sh`
- `python3 modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py`
- `python3 modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py`
- `python3 modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py`
- `bash modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh`
- `bash modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh --env-file env`
