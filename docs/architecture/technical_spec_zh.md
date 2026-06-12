# Pertussis Gene 技术总说明

## 1. 范围

本文件是面向代码和脚本的统一技术说明，目标是替代“需要在很多零散 Markdown 之间来回跳转”的阅读方式。它按真实执行边界整理当前仓库中的活跃代码文件，说明每个模块如何设计、为什么这样设计、关键结果写到哪里，以及最终哪些文件进入手稿与图表生产链。

## 2. 设计原则

- 仓库以“单一事实来源 + 分阶段责任边界”为核心。`outputs/workflow/manifest/manifest.tsv` 是跨步骤统一样本清单，而 `modules/step2_typing/outputs/bp_genotype_manifest.tsv` 是先于 PRN 结构层的标准化分型合同；每个阶段尽量把复杂逻辑放在自己的 `modules/*/bin/` 中；`workflow/bin/` 和 `workflow/lib/` 只保留跨阶段编排和公共帮助脚本。

- 活跃执行路径与审计材料显式分离。审计专用脚本、旧版手稿表格和运行快照不进入 release 分支；如需历史对照，从 git history 中恢复。

- 系统发育与 ASR 使用“双轨验证”策略。M4 使用 Gubbins + IQ-TREE，必要时附加 ClonalFrameML 与 RAxML；M5 同时产出 Fitch 与 PastML 结果，用来避免单一方法决定最终演化结论。

- 装配级判定与读段级验证分层处理。Step4 先在 assembly 层给出 `prn` 机制分类，再用 `ismapper + panISa` 进行 M6 读段验证；这样设计是因为 IS 插入和重排在组装层面常会被打断或错误折叠。

- 手稿与图表层不重复分析逻辑，而是消费已经冻结的 TSV。`manuscript/scripts/` 负责把工作流结果整理成 submission-facing 表格，`manuscript/figures/scripts/` 只负责渲染。

## 3. 端到端执行链

| 阶段 | 说明 | 主要输出 |
| --- | --- | --- |
| M0 | 建立统一 manifest、运行 readiness checks、记录版本快照。 | `outputs/workflow/manifest/`, `outputs/workflow/checkpoints/`, `outputs/workflow/versions.txt` |
| M1/M2 | 完整性检查、reads 追踪与下载计划、assembly QC、missingness 模型。 | `outputs/workflow/reads_plan/`, `outputs/workflow/assembly_qc/`, `outputs/workflow/missingness_model/` |
| M3 | 以 contig-mode Snippy 建立大树 bootstrap 输入。 | `outputs/workflow/snippy_ctg/`, `outputs/workflow/phylo/core.full.aln` |
| M4 | 缺失率预过滤、Gubbins、重组掩蔽、IQ-TREE2、可选 CFML/RAxML。 | `outputs/workflow/phylo/` |
| M5 | 参考株定根、Fitch/PastML 双轨 ASR、独立起源事件打包。 | `outputs/workflow/asr/`, `outputs/workflow/asr_sensitivity/` |
| M6 / Step4 | 对 assembly 的 `prn` 机制判定做 reads 复核并生成 hotspot 统计。 | `modules/step4_prn_validation/outputs/` |
| M7 / PH + Step6 | 构建国家-年份公共卫生主表，与基因组汇总联接后拟合生态与传播模型。 | `modules/public_health/outputs/`, `modules/step6_epi_transmission/outputs/`, `outputs/bp_step6_v2_fixed/` |
| Manuscript / Figures | 生成冻结 supplementary、figure-data extracts 和最终图件。 | `manuscript/figure_data/`, `manuscript/supplementary/`, `manuscript/figures/outputs/main/`, `manuscript/figures/outputs/extended_data/` |

## 4. 结果路径总览

- `outputs/workflow/manifest/`: cohort single source of truth, reads linkage, build report.
- `outputs/workflow/checkpoints/`: readiness JSON 报告与 consolidated readiness report.
- `outputs/workflow/reads_plan/`: read availability and download planning.
- `outputs/workflow/assembly_qc/`: assembly QC statistics and summaries.
- `outputs/workflow/missingness_model/`: interpretablity and selection-bias diagnostics.
- `outputs/workflow/snippy_ctg/`: Snippy plan, batch ledgers, completion summaries.
- `outputs/workflow/phylo/`: core alignment, pre-Gubbins filtering, recombination masks, ML trees, tree comparison.
- `outputs/workflow/asr/`: rooted tree, tip states, Fitch states, PastML states, origin events, event subtrees.
- `outputs/workflow/asr_sensitivity/`: support-threshold and composition-filtered robustness runs.
- `modules/step4_prn_validation/outputs/`: PRN mechanism calls, read validation, IS evidence, TSD, hotspot outputs.
- `modules/public_health/outputs/`: cleaned WHO/WUENIC/AMU tables and country-year master table.
- `modules/step6_epi_transmission/outputs/`: joined analytical table, model tables, cross-validation, transmission outputs.
- `outputs/bp_step4_with_ci/` and `outputs/bp_step6_v2_fixed/`: compatibility output locations still consumed by the current figure layer.
- `manuscript/figure_data/`: figure-ready frozen TSV extracts.
- `manuscript/supplementary/`: canonical supplementary tables and data files.
- `manuscript/figures/outputs/main/`, `manuscript/figures/outputs/extended_data/`: rendered main figures and Supplementary Figure source images.

## 5. 文件级说明

### 5.1 根目录编排与工作流核心

- [Snakefile](../Snakefile): 顶层依赖图，声明工作流“应该产出什么”，而不把业务细节写死在根文件里；这样设计可以把每个阶段的复杂性下放到 `rules/*.smk`。
- [config/workflow.yaml](../config/workflow.yaml): 集中管理阈值、参考序列路径和阶段输入；这样设计是为了避免参数散落在 shell、Python 和 Snakemake 三层。
- [rules/manifest.smk](../rules/manifest.smk): 把 cohort manifest、reads 可用性和 readiness checkpoint 纳入 Snakemake 依赖图；这样设计是为了让后续 QC/phylogeny 共用同一 cohort。
- [rules/reads_qc.smk](../rules/reads_qc.smk): 封装 reads 级 QC 与 MultiQC 聚合；这样设计是为了把批量 reads 质量控制作为可声明依赖。
- [rules/assembly_qc.smk](../rules/assembly_qc.smk): 把 QUAST 汇总和 missingness 诊断放进依赖图；这样设计是为了让装配质量与选择偏倚建模成为正式工件。
- [rules/snippy.smk](../rules/snippy.smk): 声明 Snippy 与 core alignment 的依赖关系；这样设计是为了让 M3 能被标准化调度。
- [rules/recomb_filter.smk](../rules/recomb_filter.smk): 把 Gubbins/重组掩蔽抽象为可重跑规则；这样设计是为了让系统发育主链显式依赖非重组 SNP。
- [rules/ml_tree.smk](../rules/ml_tree.smk): 负责 ML tree 及树比较工件；这样设计是为了把树推断与下游 ASR 解耦。
- [rules/asr.smk](../rules/asr.smk): 把定根、Fitch、PastML 和 origin packaging 放进统一依赖图；这样设计是为了把事件级结论变成可追踪工件。
- [rules/is_detection.smk](../rules/is_detection.smk): 把 Step4 的 reads 验证输出接入根工作流，但不复制 Step4 逻辑；这样设计是为了保持“stage-owned logic”。
- [rules/epi_models.smk](../rules/epi_models.smk): 把 exposure index、IPW 与 panel model 纳入依赖图；这样设计是为了让流行病学层可被重复生成。
- [workflow/bin/run_full_workflow.sh](../workflow/bin/run_full_workflow.sh): 当前最重要的根级总入口，把 foundation→M1/M2→M3→M4→M5 串成可恢复执行链，并在每一步后验证关键产物；这样设计是为了在 Snakemake 之外保留一条稳定的生产路径。
- [workflow/bin/m0_foundation.sh](../workflow/bin/m0_foundation.sh): 生成 manifest、运行 readiness checks、写版本快照；这样设计是为了把所有后续步骤的 cohort 与 readiness 决策固定下来。
- [workflow/bin/m1_m2_qc.sh](../workflow/bin/m1_m2_qc.sh): 包装完整性检查、reads 计划、assembly QC 与 missingness model；这样设计是为了把“是否可解释”与“为何不可解释”作为统一阶段。
- [workflow/bin/m3_snippy.sh](../workflow/bin/m3_snippy.sh): 用 plan→Snippy batch→snippy-core 的方式构建 bootstrap alignment；这样设计是为了支持 staged batching 和在主机资源受限时渐进扩树。
- [workflow/bin/m4_phylogeny.sh](../workflow/bin/m4_phylogeny.sh): 把缺失率过滤、Gubbins、重组掩蔽、IQ-TREE2 以及可选 CFML/RAxML 串成一条宿主友好的系统发育主链；这样设计是为了在 Conda 与 Docker 环境不稳定时仍可复跑。
- [workflow/bin/m5_asr.sh](../workflow/bin/m5_asr.sh): M5 主入口，完成参考株定根、Fitch、PastML、origin packaging 与汇总；这样设计是为了让所有 ASR 工件从同一 rooted tree 和 manifest 派生。
- [workflow/bin/m5_asr_sensitivity.sh](../workflow/bin/m5_asr_sensitivity.sh): 执行 support threshold 与 composition-filtered 两类稳健性情景；这样设计是为了把“事件数是否稳定”做成正式输出而不是口头说明。
- [workflow/lib/run_foundation_checks.py](../workflow/lib/run_foundation_checks.py): 聚合 reads availability、vaccine-variable coverage 和 validation feasibility，并统一输出 JSON 报告；这样设计是为了让 readiness 决策可自动化复跑。
- [workflow/lib/build_analysis_manifest.py](../workflow/lib/build_analysis_manifest.py): 用 Step4 机制表定义 manuscript cohort，再用 Step5/Step1 等补足 reads 和元数据；这样设计是为了让整个项目只有一个 cohort SSOT。
- [workflow/lib/aggregate_assembly_qc.py](../workflow/lib/aggregate_assembly_qc.py): 汇总 per-sample assembly QC 到统一表；这样设计是为了让 root workflow 可以消费一致的 QC 接口。
- [workflow/lib/missingness_model.py](../workflow/lib/missingness_model.py): 对 `prn` 可解释性与缺失机制建模；这样设计是为了把 selection bias 变成显式诊断层。
- [workflow/lib/filter_alignment_by_missingness.py](../workflow/lib/filter_alignment_by_missingness.py): 在 Gubbins 前过滤高缺失 tip；这样设计是为了避免极端缺失样本污染后续重组与树推断。
- [workflow/lib/mask_recombination.py](../workflow/lib/mask_recombination.py): 根据 Gubbins GFF 对 alignment 做重组位点掩蔽；这样设计是为了保留非重组 SNP 主链。
- [workflow/lib/extract_iqtree_composition_report.py](../workflow/lib/extract_iqtree_composition_report.py): 把 IQ-TREE 日志中的 composition failure 解析成 TSV；这样设计是为了让 composition-filtered sensitivity 可以自动消费。
- [workflow/lib/compare_trees.py](../workflow/lib/compare_trees.py): 比较 IQ-TREE 与 RAxML 的 tip set 和一致性；这样设计是为了把“树是否稳定”做成机器可读比较。
- [workflow/lib/root_tree_on_tip.py](../workflow/lib/root_tree_on_tip.py): 以参考株为外群定根并生成节点元数据；这样设计是为了让 Fitch 和 PastML 共用完全一致的 rooted topology。
- [workflow/lib/asr_parsimony.py](../workflow/lib/asr_parsimony.py): 实现 Fitch parsimony，并输出 tip states、node states、transitions 和 PastML 输入；这样设计是为了把审计表和下游 likelihood 输入一次性生产。
- [workflow/lib/asr_pastml.py](../workflow/lib/asr_pastml.py): 把 PastML 原始输出标准化，并与 Fitch origin 进行对照；这样设计是为了把 likelihood 结果转换为可与 parsimony 对比的结构化表。
- [workflow/lib/origin_events.py](../workflow/lib/origin_events.py): 扫描 intact→disrupted 转变并打包 descendant tip 子集；这样设计是为了把“独立起源”从概念变成事件级工件。
- [workflow/lib/prune_tree_by_tips.py](../workflow/lib/prune_tree_by_tips.py): 按 tip 列表修剪 Newick 树；这样设计是为了不依赖 `ete3` 也能完成 composition-filtered sensitivity。
- [workflow/lib/build_ap_exposure_index.py](../workflow/lib/build_ap_exposure_index.py): 把疫苗制度变量组合成 aP exposure index；这样设计是为了把“免疫暴露”从原始表提升为模型可消费的标准化输入。
- [workflow/lib/ipw_prevalence.py](../workflow/lib/ipw_prevalence.py): 对 `prn` disrupted 频率做 IPW 校正；这样设计是为了把 interpretability bias 显式纳入生态分析。
- [workflow/lib/panel_model.py](../workflow/lib/panel_model.py): 执行国家-年份面板模型；这样设计是为了把 exposure index 与 `prn` disrupted 频率联系起来。
- [workflow/lib/its_feasibility.py](../workflow/lib/its_feasibility.py): 检查 interrupted time series 是否具备最基本可行性；这样设计是为了在样本和年份不足时及时阻断伪因果叙述。
- [manuscript/scripts/freeze/generate_supplementary_table_1.py](../manuscript/scripts/freeze/generate_supplementary_table_1.py): 导出手稿 Supplementary Table 1；这样设计是为了避免手工整理 cohort metadata。
- [workflow/bin/find-big-git-objects.sh](../workflow/bin/find-big-git-objects.sh): 仓库维护脚本，用于定位超大 Git 对象；这样设计是为了防止分析仓库再次被大文件污染。
- [validate_figure1.py](../validate_figure1.py): 一次性校验旧版与修正版 `R_e` 数据质量差异；这样设计是为了在迁移到 `outputs/bp_step6_v2_fixed/` 后快速确认 Figure 1 输入已修正。

### 5.2 Manifest / Raw Reads

- [modules/step1_ingest/bin/core/01_fetch_ncbi_report.py](../modules/step1_ingest/bin/core/01_fetch_ncbi_report.py): 抓取 NCBI 原始 metadata 报告；设计上把“下载原始元数据”与“清洗/聚合”拆开，便于重抓取。
- [modules/step1_ingest/bin/core/02_export_ncbi_tsv.py](../modules/step1_ingest/bin/core/02_export_ncbi_tsv.py): 把 NCBI JSON/报告展平成 TSV；这样设计是为了给后续清洗脚本提供稳定表结构。
- [modules/step1_ingest/bin/core/03_clean_metadata_aggregate.py](../modules/step1_ingest/bin/core/03_clean_metadata_aggregate.py): 标准化国家、日期和来源字段并聚合；这样设计是为了把 source heterogeneity 在进入 cohort 前处理掉。
- [modules/step1_ingest/bin/core/04_download_ncbi_genomes.py](../modules/step1_ingest/bin/core/04_download_ncbi_genomes.py): 下载 assembly FASTA；这样设计是为了把 metadata 与实体基因组同步到本地/NAS。
- [modules/step1_ingest/bin/core/05_verify_genome_extract.py](../modules/step1_ingest/bin/core/05_verify_genome_extract.py): 核对下载完整性与提取结果；这样设计是为了在进入 QC 之前先拦截下载损坏。
- [modules/step1_ingest/bin/manifest/06_build_public_manifest.py](../modules/step1_ingest/bin/manifest/06_build_public_manifest.py): 建立未去重的公开样本清单；这样设计是为了先保留最大覆盖，再在后续显式决策去重。
- [modules/step1_ingest/bin/manifest/07_recover_raw_reads.py](../modules/step1_ingest/bin/manifest/07_recover_raw_reads.py): 从 SRA/ENA 追溯 reads 关联；这样设计是为了把 assembly-only 与 read-backed 路径在 manifest 层连起来。
- [modules/step1_ingest/bin/manifest/08_resolve_duplicates.py](../modules/step1_ingest/bin/manifest/08_resolve_duplicates.py): 解决重复 accession/样本记录；这样设计是为了给下游模型和 phylogeny 提供唯一 canonical sample。
- [modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py](../modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py): 构建分析 cohort，尤其是后续 country-year 层使用的 cohort C；这样设计是为了把不同分析边界显式编码。
- [modules/step1_ingest/bin/raw_reads/10_build_download_plan.py](../modules/step1_ingest/bin/raw_reads/10_build_download_plan.py): 生成原始 reads 下载计划；这样设计是为了让重度 IO 操作可以批处理。
- [modules/step1_ingest/bin/raw_reads/11_split_download_plan.py](../modules/step1_ingest/bin/raw_reads/11_split_download_plan.py): 把大下载计划拆成 shard；这样设计是为了支持分布式或断点续跑。
- [modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh](../modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh): 执行单个 shard 的下载/组装链；这样设计是为了让大规模组装任务能按分片执行。
- [modules/step1_ingest/bin/raw_reads/13_preflight_env.sh](../modules/step1_ingest/bin/raw_reads/13_preflight_env.sh): 环境预检；这样设计是为了在正式下载和组装前尽早发现工具缺失。
- [modules/step1_ingest/bin/raw_reads/14_setup_env.sh](../modules/step1_ingest/bin/raw_reads/14_setup_env.sh): 准备运行环境；这样设计是为了标准化多机或多目录的执行前置条件。
- [modules/step1_ingest/bin/raw_reads/15_fetch_taxon_read_run_catalog.py](../modules/step1_ingest/bin/raw_reads/15_fetch_taxon_read_run_catalog.py): 抓取 taxon 级 reads/run catalog；这样设计是为了从 accession 层上游补足 run 候选。
- [modules/step1_ingest/bin/raw_reads/16_build_external_gapfill.py](../modules/step1_ingest/bin/raw_reads/16_build_external_gapfill.py): 为缺口样本构建外部补齐计划；这样设计是为了在现有 assembly/reads 不足时保留补救路径。
- [modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py](../modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py): 收集分片组装结果；这样设计是为了把分布式产出重新并回统一目录。
- [modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py](../modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py): 对新组装结果做 QC；这样设计是为了在汇入主数据库前先筛掉明显失败样本。
- [modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py](../modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py): 合并通过 QC 的组装结果；这样设计是为了把下载/组装侧的结果与主 assembly 库同步。
- [modules/step1_ingest/bin/raw_reads/20_build_genome_database.py](../modules/step1_ingest/bin/raw_reads/20_build_genome_database.py): 构建 `pertussis_data/bp_genomes_qc/assemblies/` 统一 genome database；这样设计是为了让 M3/M4/M5 都只面对一个 assembly 根目录。
- [modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh](../modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh): 下载缺失 assembly；这样设计是为了把“数据库补齐”做成独立操作。
- [modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py](../modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py): 检查 genome database 完整性；这样设计是为了让 M1 先回答“数据是否齐全”。
- [modules/step1_ingest/bin/raw_reads/23_retry_missing_assemblies.sh](../modules/step1_ingest/bin/raw_reads/23_retry_missing_assemblies.sh): 重试失败下载；这样设计是为了把网络失败和永久缺失分开。
- [modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py](../modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py): 根据 manifest 追踪样本是否有原始 reads；这样设计是为了支撑 readiness checkpoint 和后续 read validation。
- [modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py](../modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py): 生成正式 reads 下载计划；这样设计是为了从“可追溯”过渡到“可执行下载”。
- [modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py](../modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py): 执行 assembly QC；这样设计是为了让装配质量成为下游 Snippy 与 missingness model 的正式输入。
- [modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh](../modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh): 运行批量 contig-mode Snippy；这样设计是为了支持大规模 staged batching。
- [modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py](../modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py): 构建 Snippy 计划表；这样设计是为了把 include/exclude 决策显式表格化。
- [modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh](../modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh): 聚合完成的 Snippy 目录并生成 core alignment；这样设计是为了让 M3 既支持 batch 也支持 cumulative rebuild。

### 5.3 Step2 装配表征

- [modules/step2_typing/bin/step2_01_qc_filter.py](../modules/step2_typing/bin/step2_01_qc_filter.py): 对 assembly 做长度、contig 数和 N50 过滤；这样设计是为了把低质量样本挡在 marker 扫描之前。
- [modules/step2_typing/bin/step2_02_index_genomes.py](../modules/step2_typing/bin/step2_02_index_genomes.py): 索引待处理基因组；这样设计是为了给批量 BLAST/MLST 提供统一输入列表。
- [modules/step2_typing/bin/step2_03_run_mlst.py](../modules/step2_typing/bin/step2_03_run_mlst.py): 运行 MLST；这样设计是为了把序列类型作为下游 lineage/summary 的标准字段。
- [modules/step2_typing/bin/step2_04_marker_scan_blast.py](../modules/step2_typing/bin/step2_04_marker_scan_blast.py): 用 BLAST 扫描关键 marker；这样设计是为了从 assembly 侧统一获得 `prn`、23S 等位点信号。
- [modules/step2_typing/bin/step2_05_merge_qc_tables.py](../modules/step2_typing/bin/step2_05_merge_qc_tables.py): 合并 QC 中间表；这样设计是为了让 MLST 与 marker 层面对统一 QC 输入。
- [modules/step2_typing/bin/step2_06_merge_mlst.py](../modules/step2_typing/bin/step2_06_merge_mlst.py): 汇总 MLST 结果；这样设计是为了把 per-sample `mlst` 输出整合成主表。
- [modules/step2_typing/bin/step2_07_mlst_summaries.py](../modules/step2_typing/bin/step2_07_mlst_summaries.py): 生成 MLST 汇总统计；这样设计是为了把 ST 分布从明细表提升为 cohort 描述。
- [modules/step2_typing/bin/step2_08_extract_marker_alleles.py](../modules/step2_typing/bin/step2_08_extract_marker_alleles.py): 提取 marker allele 序列与哈希；这样设计是为了把“已知/未知 allele”都放进统一比较框架。
- [modules/step2_typing/bin/step2_09_call_23s_a2047g.py](../modules/step2_typing/bin/step2_09_call_23s_a2047g.py): 专门识别 23S A2047G；这样设计是为了把耐大环内酯标志位点做成显式字段。
- [modules/step2_typing/bin/step2_10_merge_markers.py](../modules/step2_typing/bin/step2_10_merge_markers.py): 合并 marker hits 到统一表；这样设计是为了让下游步骤只消费一张 marker 主表。
- [modules/step2_typing/bin/step2_11_marker_summaries.py](../modules/step2_typing/bin/step2_11_marker_summaries.py): 生成 marker 汇总；这样设计是为了给 Step3/Step6 提供更轻量的摘要接口。
- [modules/step2_typing/bin/step2_12_build_marker_references.py](../modules/step2_typing/bin/step2_12_build_marker_references.py): 生成 marker reference 资产；这样设计是为了保证 BLAST/allele extraction 使用一致参考。
- [modules/step2_typing/bin/step2_13_joint_summaries.py](../modules/step2_typing/bin/step2_13_joint_summaries.py): 把 QC、MLST 与 marker 结果拼接成 Step2 主表；这样设计是为了让 Step4/Step6 只读取一张综合表。
- [modules/step2_typing/bin/step2_14_harmonize_typing.py](../modules/step2_typing/bin/step2_14_harmonize_typing.py): 把原始 marker hash 统一到标准化 `ptxP`/`fim3`/`fhaB2400_5550`/`23S` 命名，并联接冻结的 typing profile registry；这样设计是为了让统一 manifest 和手稿层消费“标准化分型层 + PRN 结构层”的双层合同。输出：`modules/step2_typing/outputs/bp_genotype_manifest.tsv`。状态：active。

### 5.4 Step3 预备层

- [modules/step3_prn_scan/bin/step3_01_extra_summaries.py](../modules/step3_prn_scan/bin/step3_01_extra_summaries.py): 生成额外 cohort 摘要；这样设计是为了在大树与机制分析前先有描述性统计。
- [modules/step3_prn_scan/bin/step3_10_prepare_phylogeny_manifest.py](../modules/step3_prn_scan/bin/step3_10_prepare_phylogeny_manifest.py): 构建用于大树的 phylogeny manifest；这样设计是为了控制大树抽样和平衡。
- [modules/step3_prn_scan/bin/step3_20_prn_disruption_scan.py](../modules/step3_prn_scan/bin/step3_20_prn_disruption_scan.py): 初步扫描 `prn` disruption；这样设计是为了从 assembly 层快速识别候选异常样本。
- [modules/step3_prn_scan/bin/step3_21_prn_disruption_summaries.py](../modules/step3_prn_scan/bin/step3_21_prn_disruption_summaries.py): 汇总 Step3 disruption 结果；这样设计是为了在进入 Step4 前先得到 cohort-level 轮廓。
- [modules/step3_prn_scan/bin/step3_30_prn_trends_tables.py](../modules/step3_prn_scan/bin/step3_30_prn_trends_tables.py): 构建 `prn` 趋势表；这样设计是为了给时序描述和后续公共卫生整合提供早期接口。
- [modules/step3_prn_scan/bin/step3_40_phylo_annotations.py](../modules/step3_prn_scan/bin/step3_40_phylo_annotations.py): 为 phylogeny manifest 添加注释；这样设计是为了把机制与谱系标签提前对齐。
- [modules/step3_prn_scan/bin/step3_50_prn_breakpoint_evidence.py](../modules/step3_prn_scan/bin/step3_50_prn_breakpoint_evidence.py): 提取断点证据；这样设计是为了在 Step4 reads 验证前先从 assembly 层收集重排线索。
- [modules/step3_prn_scan/bin/step3_51_prn_breakpoint_summaries.py](../modules/step3_prn_scan/bin/step3_51_prn_breakpoint_summaries.py): 汇总断点证据；这样设计是为了把个别样本证据转成 cohort 级摘要。
- [modules/step3_prn_scan/bin/step3_52_extract_prn_gap_sequences.py](../modules/step3_prn_scan/bin/step3_52_extract_prn_gap_sequences.py): 提取 `prn` gap 区序列；这样设计是为了支持 Step4 IS 扫描和 gap 解释。
- [modules/step3_prn_scan/bin/step3_60_results_digest.py](../modules/step3_prn_scan/bin/step3_60_results_digest.py): 打包 Step3 结果摘要；这样设计是为了给手动检查和后续步骤提供轻量汇总。

### 5.5 Step4 机制解析与读段验证

- [modules/step4_prn_validation/bin/step4_00_distributed_raw_reads_lib.sh](../modules/step4_prn_validation/bin/step4_00_distributed_raw_reads_lib.sh): 分布式 raw-read 作业的公共 shell 库；这样设计是为了把多机执行的共用逻辑抽离。
- [modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh](../modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh): 启动分布式 reads 作业；这样设计是为了让大批量 Step4 任务可以分片派发。
- [modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh](../modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh): 收集分布式作业状态；这样设计是为了统一查看长时间运行任务。
- [modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh](../modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh): 同步分布式输出；这样设计是为了把多节点产物收回主仓库。
- [modules/step4_prn_validation/bin/step4_00_rebalance_distributed_shards.py](../modules/step4_prn_validation/bin/step4_00_rebalance_distributed_shards.py): 重平衡分片；这样设计是为了避免不同节点负载过于失衡。
- [modules/step4_prn_validation/bin/step4_01_build_is_reference.py](../modules/step4_prn_validation/bin/step4_01_build_is_reference.py): 构建 IS reference 资产；这样设计是为了让 assembly 扫描和 reads 验证共享同一 IS 参考。
- [modules/step4_prn_validation/bin/step4_02_scan_prn_mechanisms.py](../modules/step4_prn_validation/bin/step4_02_scan_prn_mechanisms.py): 在 assembly 层扫描 `prn` 机制；这样设计是为了先大规模获得机制候选，再用 reads 做重点复核。
- [modules/step4_prn_validation/bin/step4_02b_summarize_is_hits.py](../modules/step4_prn_validation/bin/step4_02b_summarize_is_hits.py): 汇总 IS hits；这样设计是为了把局部比对证据提升为机制层摘要。
- [modules/step4_prn_validation/bin/step4_02c_score_prn_calls.py](../modules/step4_prn_validation/bin/step4_02c_score_prn_calls.py): 为 `prn` 判定打置信度分层；这样设计是为了区分 assembly-high/moderate/low 与不足证据。
- [modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables.py](../modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables.py): 旧版 summary table 生成器；保留用于兼容旧输出。
- [modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables_v2.py](../modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables_v2.py): 当前主用 summary 生成器，加入更稳定的 manuscript-facing 汇总逻辑与置信区间；这样设计是为了支持当前图表层。
- [modules/step4_prn_validation/bin/step4_03a_build_validation_subset.py](../modules/step4_prn_validation/bin/step4_03a_build_validation_subset.py): 构建 reads 验证子集；这样设计是为了把稀缺 reads 资源优先用于最有信息量的样本。
- [modules/step4_prn_validation/bin/step4_03b_assess_validation_feasibility.py](../modules/step4_prn_validation/bin/step4_03b_assess_validation_feasibility.py): validation feasibility 检查；这样设计是为了在启动 M6 前先确认是否有足够 reads 候选。
- [modules/step4_prn_validation/bin/step4_03c_prepare_is_reference.py](../modules/step4_prn_validation/bin/step4_03c_prepare_is_reference.py): 为 reads 验证准备专门 IS reference；这样设计是为了让 `ismapper/panISa` 输入保持一致。
- [modules/step4_prn_validation/bin/step4_03d_build_read_validation_batch.py](../modules/step4_prn_validation/bin/step4_03d_build_read_validation_batch.py): 把验证子集转成批处理计划；这样设计是为了让下载、比对和验证在一个 batch contract 上运行。
- [modules/step4_prn_validation/bin/step4_03e_run_is_read_validation.sh](../modules/step4_prn_validation/bin/step4_03e_run_is_read_validation.sh): M6 主 shell 入口，调度 `ismapper`、`panISa` 与下游解析；这样设计是为了把长链路 reads 验证包成一个可操作的批处理入口。
- [modules/step4_prn_validation/bin/step4_03_validate_prn_with_reads.py](../modules/step4_prn_validation/bin/step4_03_validate_prn_with_reads.py): 读段验证核心解析器，把 `ismapper` 与 `panISa` 证据合并为 `supported / supported_concordant / unresolved / no_signal` 等状态；这样设计是为了把多工具证据压缩成一个手稿可用的验证状态层。
- [modules/step4_prn_validation/bin/step4_03f_hotspot_test.py](../modules/step4_prn_validation/bin/step4_03f_hotspot_test.py): 检验插入热点；这样设计是为了把“某个位点是否聚集”转化为统计工件。
- [modules/step4_prn_validation/bin/step4_04_summarize_prn_validation.py](../modules/step4_prn_validation/bin/step4_04_summarize_prn_validation.py): 汇总 reads 验证结果；这样设计是为了让 Step4 输出既有明细也有 cohort 摘要。
- [modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py](../modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py): 构建当前 M1/M2 wrapper 使用的 missingness model；这样设计是为了把缺失机制建模放回已验证的根编排路径。输出：`outputs/workflow/missingness_model/`。状态：active dependency。

### 5.6 Step5 旧平衡树层

- [modules/step5_phylogeny_asr/bin/step5_01_build_phylogeny_manifest.py](../modules/step5_phylogeny_asr/bin/step5_01_build_phylogeny_manifest.py): 构建 Step5 旧版 phylogeny manifest；这样设计是为了在 workflow-native M3/M4/M5 之前先支持大树抽样。
- [modules/step5_phylogeny_asr/bin/step5_02_build_global_phylogeny.sh](../modules/step5_phylogeny_asr/bin/step5_02_build_global_phylogeny.sh): 旧版全局树构建入口；保留供历史结果追溯。
- [modules/step5_phylogeny_asr/bin/step5_03_reconstruct_prn_states.py](../modules/step5_phylogeny_asr/bin/step5_03_reconstruct_prn_states.py): 旧版 ASR 状态重建；这样设计是为了在 workflow-native ASR 之前提供平衡树上的起源计数。
- [modules/step5_phylogeny_asr/bin/step5_04_count_independent_origins.py](../modules/step5_phylogeny_asr/bin/step5_04_count_independent_origins.py): 旧版独立起源计数；这样设计是为了把平衡树层的事件数做成表格。
- [modules/step5_phylogeny_asr/bin/step5_05_summarize_clades.py](../modules/step5_phylogeny_asr/bin/step5_05_summarize_clades.py): 汇总 clade 层结果；这样设计是为了给旧版系统发育分析提供摘要层。
- [modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py](../modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py): 当前 M1/M2 仍直接调用的 missingness model 脚本；这样设计是因为它在根编排中已被验证。
- [modules/step5_phylogeny_asr/bin/step5_07_asr_sensitivity.py](../modules/step5_phylogeny_asr/bin/step5_07_asr_sensitivity.py): 旧版 ASR sensitivity；这样设计是为了在 balanced-tree 时代提供稳健性情景。

### 5.7 公共卫生数据层

- [modules/public_health/bin/ph_utils.py](../modules/public_health/bin/ph_utils.py): 公共卫生层共享工具，负责文本标准化、Excel 读取和 TSV 写出；这样设计是为了避免每个清洗脚本重复处理国家名和 freeze date。
- [modules/public_health/bin/ph_01_build_source_inventory.py](../modules/public_health/bin/ph_01_build_source_inventory.py): 建立源数据清单；这样设计是为了把所有输入来源版本化。
- [modules/public_health/bin/ph_02_normalize_country_names.py](../modules/public_health/bin/ph_02_normalize_country_names.py): 标准化国家名称与 ISO3；这样设计是为了让 WHO/WUENIC/AMU/基因组层能按国家-年份联接。
- [modules/public_health/bin/ph_03_clean_who_cases.py](../modules/public_health/bin/ph_03_clean_who_cases.py): 清洗 WHO 病例数据；这样设计是为了把病例表转成可与国家-年份面板直接联接的形式。
- [modules/public_health/bin/ph_04_clean_wuenic.py](../modules/public_health/bin/ph_04_clean_wuenic.py): 清洗 WUENIC 免疫覆盖率数据；这样设计是为了把 DTP3/booster 等字段统一进时间序列框架。
- [modules/public_health/bin/ph_05_clean_vaccine_programs.py](../modules/public_health/bin/ph_05_clean_vaccine_programs.py): 整理疫苗项目元数据，包括 acellular/whole-cell 与 `prn_in_vaccine`；这样设计是为了把“暴露定义”做成显式字段。
- [modules/public_health/bin/ph_06_clean_glass_amu.py](../modules/public_health/bin/ph_06_clean_glass_amu.py): 清洗 GLASS 抗菌药使用数据；这样设计是为了与 ESAC-Net 形成互补来源。
- [modules/public_health/bin/ph_07_clean_esacnet_amu.py](../modules/public_health/bin/ph_07_clean_esacnet_amu.py): 清洗 ESAC-Net AMU 数据；这样设计是为了在 GLASS 缺失时补足 macrolide/total antibiotic 指标。
- [modules/public_health/bin/ph_08_build_country_year_master.py](../modules/public_health/bin/ph_08_build_country_year_master.py): 把病例、免疫项目和 AMU 合并成统一国家-年份主表；这样设计是为了给 Step6 提供唯一公共卫生输入层。
- [modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py](../modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py): 检查疫苗变量是否足够支持后续生态模型，并默认读取更丰富的 formulation curation 文件；这样设计是为了在暴露变量不足时提前降级叙述。

### 5.8 Step6 生态与传播建模

- [modules/step6_epi_transmission/bin/step6_01_build_country_year_genomic_summaries.py](../modules/step6_epi_transmission/bin/step6_01_build_country_year_genomic_summaries.py): 把样本级 `prn` 机制、reads 支持和 A2047G 汇总到国家-年份；这样设计是为了把 genome-level 结果投影到公共卫生分析层。
- [modules/step6_epi_transmission/bin/step6_02_join_public_health.py](../modules/step6_epi_transmission/bin/step6_02_join_public_health.py): 将基因组汇总与 `ph_country_year_master.tsv` 合并；这样设计是为了形成模型统一输入表。
- [modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py](../modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py): 拟合主生态模型；这样设计是为了给手稿主效应图提供标准输出，并在样本稀疏时显式声明未做随机效应。
- [modules/step6_epi_transmission/bin/step6_04_run_sensitivity_models.py](../modules/step6_epi_transmission/bin/step6_04_run_sensitivity_models.py): 运行生态模型敏感性分析；这样设计是为了比较不同 covariate/过滤策略下的稳定性。
- [modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py](../modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py): 执行 AMU 探索性敏感性分析；这样设计是为了把 AMU 信号保持在“探索性”而非主结论层。
- [modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers.py](../modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers.py): 第一代 `R_e` 估计实现；保留用于追溯传播动力学开发过程。
- [modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers_v2.py](../modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers_v2.py): 修正版 `R_e` 估计实现；这样设计是为了清洗异常值并给图件提供质量标记。
- [modules/step6_epi_transmission/bin/step6_06_run_re_estimation.sh](../modules/step6_epi_transmission/bin/step6_06_run_re_estimation.sh): 包装 `R_e` 估计运行；这样设计是为了把输入输出目录和运行参数标准化。
- [modules/step6_epi_transmission/bin/step6_07_fit_transmission_models.py](../modules/step6_epi_transmission/bin/step6_07_fit_transmission_models.py): 拟合传播模型；这样设计是为了把 `R_e` 轨迹与国家层 covariates 联系起来。
- [modules/step6_epi_transmission/bin/step6_07_run_transmission_models.sh](../modules/step6_epi_transmission/bin/step6_07_run_transmission_models.sh): 传播模型 shell 包装器；这样设计是为了让 `R_e` 输入和模型输出有固定接口。
- [modules/step6_epi_transmission/bin/step6_08_cross_validation.py](../modules/step6_epi_transmission/bin/step6_08_cross_validation.py): 做模型交叉验证；这样设计是为了让探索性模型仍有基本泛化诊断。
- [modules/step6_epi_transmission/bin/step6_09_mixed_effects_models.py](../modules/step6_epi_transmission/bin/step6_09_mixed_effects_models.py): 尝试 mixed-effects 版本；这样设计是为了在数据允许时提供更强结构化模型。

### 5.9 手稿冻结层

- [manuscript/scripts/freeze/extract_key_statistics.py](../manuscript/scripts/freeze/extract_key_statistics.py): 抽取手稿关键数字并写入 JSON/TXT；这样设计是为了避免摘要和结果部分手工抄写。
- [manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py](../manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py): 把 Step4、Step6、ASR 等结果切成 figure-ready TSV，并同时生成数据字典；这样设计是为了让图件层只消费冻结抽取结果。
- [manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py](../manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py): 把 workflow-native ASR 输出复制到 supplementary 和 figure-data 目录，并把 composition-pruned ASR 框架升为主分析；这样设计是为了把主工作流结果显式升格为 submission-facing 工件。
- [manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py](../manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py): 把 composition-pruned rooted ML tree 转成 Figure 3 可直接绘制的 segments/nodes 表；这样设计是为了让 Figure 3 的树图完全脱离 R 端的在线树解析。
- [manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py](../manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py): 诊断总调度脚本，统一触发 ASR、验证证据和上下文诊断输出。
- [manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py](../manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py): 生成合并后的 ASR、validation 和 context 证据层；输出 `Supplementary_Table_13` 到 `Supplementary_Table_19`、`Supplementary_Table_21`、`Supplementary_Table_23` 和 `Supplementary_Table_24`，以及对应的 Figure 3 和 validation sidecars。
- [manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py](../manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py): 生成 submission evidence summary 以及 lineage/origin collapse 表；输出 `Supplementary_Table_25` 到 `Supplementary_Table_30` 和 consolidated validation matrix。
- [manuscript/scripts/review/ms_15_build_selected_country_review_report.py](../manuscript/scripts/review/ms_15_build_selected_country_review_report.py): 生成 selected-country 的 program-history manifest、分层表和综合证据网格；输出 `manuscript/figure_data/selected_country/*`。
- [manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py](../manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py): 生成四个 reviewer-facing sidecar 表，覆盖年份敏感性、architecture turnover、origin burden 连接和 PRN-specificity 负对照；输出 `manuscript/figure_data/selected_country/*` 以及 `manuscript/supplementary/Supplementary_Table_35_*.tsv` 到 `Supplementary_Table_38_*.tsv`。
- [manuscript/scripts/source_data/ms_14_build_source_data_manifest.py](../manuscript/scripts/source_data/ms_14_build_source_data_manifest.py): 生成当前 Figure 1-5 与 Extended Data Fig. 1-10 合同对应的 source-data 清单；输出 `manuscript/submission_data/source_data/final_source_data_manifest.tsv`。
- [manuscript/scripts/source_data/ms_17_build_source_data_workbook.py](../manuscript/scripts/source_data/ms_17_build_source_data_workbook.py): 基于 source-data manifest 生成 panel/file 清单与逐图 Excel workbook；输出 `manuscript/submission_data/source_data/`。

### 5.10 图件渲染层

- [manuscript/figures/bin/render_main.R](../manuscript/figures/bin/render_main.R): 主图总入口，顺序调用 Figure 1-5；这样设计是为了让主文图件一键复跑。
- [manuscript/figures/bin/render_extended_data.R](../manuscript/figures/bin/render_extended_data.R): Extended Data 总入口，顺序调用 Extended Data Fig. 1-10；这样设计是为了把主文和 ED 的生成边界彻底分开。
- [manuscript/figures/scripts/main/fig01_public_genome_atlas.R](../manuscript/figures/scripts/main/fig01_public_genome_atlas.R): 渲染 public-genome atlas、采样框架和可判读 *prn* 边界；这样设计是为了先固定数据宇宙与 locus recoverability。
- [manuscript/figures/scripts/main/fig02_prn_structural_solution_space.R](../manuscript/figures/scripts/main/fig02_prn_structural_solution_space.R): 渲染成比例的 *prn* locus atlas 与结构事件复用面板；这样设计是为了把结构约束前置为机制主线。
- [manuscript/figures/scripts/main/fig03_repeated_origin_phylogeny.R](../manuscript/figures/scripts/main/fig03_repeated_origin_phylogeny.R): 渲染带外圈上下文轨道的 fan tree repeated-origin 分析和 ASR sensitivity 面板；这样设计是为了把 phylogeny 与验证上下文合成在一个主图。
- [manuscript/figures/scripts/main/fig04_country_programme_amplification.R](../manuscript/figures/scripts/main/fig04_country_programme_amplification.R): 渲染国家疫苗项目 amplification 小多图和 selected-country trajectory contrasts；这样设计是为了把疫苗项目史作为放大环境而非旧的图件主轴。
- [manuscript/figures/scripts/main/fig05_validation_synthesis.R](../manuscript/figures/scripts/main/fig05_validation_synthesis.R): 渲染 PRN specificity、support-only exposure sensitivity、missingness contrasts 和 identifiability synthesis grid；这样设计是为了用验证和可识别性边界收束主文。
- [manuscript/figures/scripts/lib/data_utils.R](../manuscript/figures/scripts/lib/data_utils.R): 图件层数据读取与校验库；现在会动态解析仓库根目录并统一读取冻结后的 manuscript-facing TSV，这样设计是为了把图件脚本从脆弱的工作目录假设里解耦出来。
- [manuscript/figures/scripts/lib/theme_nature.R](../manuscript/figures/scripts/lib/theme_nature.R): 统一主图主题、色板与保存函数；这样设计是为了让所有主图共享一致视觉规范。
- [manuscript/figures/scripts/extended_data/ed09_ecology_sidecar_impl.R](../manuscript/figures/scripts/extended_data/ed09_ecology_sidecar_impl.R): 绘制 support-only ecology/programme 诊断并映射为 Supplementary Figure 9 的源图资产；这样设计是为了把生态层明确降级为辅助桥接而非主结论。

### 5.11 测试覆盖

- [tests/test_m4_utils.py](../tests/test_m4_utils.py): 测试 `mask_recombination.py` 与 `compare_trees.py` 的 CLI 行为；这样设计是为了保护 M4 两个最关键纯 Python 辅助脚本。
- [tests/test_asr_pastml_parser.py](../tests/test_asr_pastml_parser.py): 测试 `asr_pastml.py` 是否能正确区分 strict 与 compatible origins；这样设计是为了保证 PastML 轨的解释稳定。
- [tests/test_m5_asr_stage1.py](../tests/test_m5_asr_stage1.py): 测试 `asr_parsimony.py` 与 `origin_events.py` 是否能从小树上恢复单个起源事件；这样设计是为了守住 M5 核心逻辑。
- [tests/test_step4_read_validation.py](../tests/test_step4_read_validation.py): 测试 `step4_03_validate_prn_with_reads.py` 与 `step4_03f_hotspot_test.py`；这样设计是为了保护 Step4 最关键的 reads 证据整合逻辑。

## 6. 推荐阅读顺序

如果要真正理解当前活跃生产链，建议按下面顺序阅读代码：

1. `workflow/bin/run_full_workflow.sh`
2. `workflow/lib/build_analysis_manifest.py` + `workflow/lib/run_foundation_checks.py`
3. `modules/step1_ingest/bin/raw_reads/24-29_*`
4. `workflow/bin/m4_phylogeny.sh` + `workflow/bin/m5_asr.sh`
5. `modules/step4_prn_validation/bin/step4_02*` 和 `step4_03*`
6. `modules/public_health/bin/ph_03-08_*`
7. `modules/step6_epi_transmission/bin/step6_01-07_*`
8. `manuscript/scripts/freeze/ms_01-03_*` + `manuscript/figures/bin/render_main.R`

## 7. 当前技术债

- 图件层目前仍依赖 `manuscript/figure_data/` 下较宽的冻结 TSV 合同；后续如果进一步收敛 staging manifest 和公共路径辅助函数，下一次渲染层迁移会更轻。

- Step5 旧 balanced-tree 逻辑仍部分参与 manifest 补充，因此还不能简单删除整个 `modules/step5_phylogeny_asr/`。

- 根级 Snakemake scaffold 与 shell 主路径并存，虽然边界已经明确，但仍属于“双入口”架构。

- 分布式 Step4 脚本和历史 Step3/Step5 结果仍保留，原因是它们对审计和再现旧分析有价值。
