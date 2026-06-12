# Pertussis project environments

当前仓库不再依赖脚本里写死的 Conda 绝对路径。

统一运行时配置在：

- `config/runtime/runtime_envs.env`: 仓库共享默认值
- `config/runtime/runtime_envs.local.env`: 本机覆盖，不进 git
- `workflow/bin/bootstrap_runtime_envs.sh`: 检查三套环境路径是否可用
- `workflow/bin/run_with_project_env.sh`: 按脚本头部声明自动选择正确环境

推荐的最小复现步骤是：

```bash
cp config/runtime/runtime_envs.local.env.example config/runtime/runtime_envs.local.env
$EDITOR config/runtime/runtime_envs.local.env
bash workflow/bin/bootstrap_runtime_envs.sh --check
```

脚本执行建议统一走：

```bash
bash workflow/bin/run_with_project_env.sh --script workflow/bin/m0_foundation.sh
bash workflow/bin/run_with_project_env.sh --script workflow/bin/m4_phylogeny.sh -- --threads 8 --iq-threads 16
```

这三个环境定义文件没有删除，只是从仓库根目录迁到了 `config/env/`：

- `config/env/environment_tool.yml` -> `pertussis-bio-tools`
- `config/env/environment_python.yml` -> `pertussis-prn-global-bio`
- `config/env/environment_r.yml` -> `pertussis-prn-global-r`

其中：

- `pertussis-bio-tools` 使用 `Python 3.10`，放绝大多数常规工具和 Python 分析依赖。
- `pertussis-prn-global-bio` 使用 `Python 3.9`，放系统发育和 Step4 读段验证这一组需要兼容性隔离的包：
  `gubbins`、`iqtree`、`raxml-ng`、`clonalframeml`、`pastml`、`ete3`、`ismapper`、`panisa`
- `pertussis-prn-global-r` 是新增的独立 R 环境。

这样拆的目的很明确：

- 需要 `Python 3.9` 的那组系统发育/兼容性包单独隔离
- 其余绝大多数包统一放到 `Python 3.10`

## Create or update with mamba

先加载运行时配置：

```bash
source config/runtime/runtime_envs.env
if [[ -f config/runtime/runtime_envs.local.env ]]; then
  source config/runtime/runtime_envs.local.env
fi
```

首次创建：

```bash
"${PROJECT_CONDA_EXE}" env create -p "${PROJECT_ENV_BIO_TOOLS_PREFIX}" -f config/env/environment_tool.yml
"${PROJECT_CONDA_EXE}" env create -p "${PROJECT_ENV_PHYLO_PREFIX}" -f config/env/environment_python.yml
"${PROJECT_CONDA_EXE}" env create -p "${PROJECT_ENV_R_PREFIX}" -f config/env/environment_r.yml
```

如果环境已经存在，推荐直接更新并清理无用依赖：

```bash
"${PROJECT_CONDA_EXE}" env update -p "${PROJECT_ENV_BIO_TOOLS_PREFIX}" -f config/env/environment_tool.yml --prune
"${PROJECT_CONDA_EXE}" env update -p "${PROJECT_ENV_PHYLO_PREFIX}" -f config/env/environment_python.yml --prune
"${PROJECT_CONDA_EXE}" env update -p "${PROJECT_ENV_R_PREFIX}" -f config/env/environment_r.yml --prune
```

如果当前 shell 里还没有加载 Conda shell hook，可以执行：

```bash
source "${PROJECT_CONDA_ROOT}/etc/profile.d/conda.sh"
```

或者完全不激活环境，直接用 launcher：

```bash
bash workflow/bin/run_with_project_env.sh --script workflow/lib/build_analysis_manifest.py -- --help
```

当前默认 channel 已切回官方 `conda-forge` / `bioconda`。原因是这台机器上的 TUNA 镜像在 `pastml` / `ete3` / `numpy` 这一组求解时出现了不同步现象，导致 `config/env/environment_python.yml` 会误报不可解。

## Activate

```bash
conda activate "${PROJECT_ENV_BIO_TOOLS_PREFIX}"
conda activate "${PROJECT_ENV_PHYLO_PREFIX}"
conda activate "${PROJECT_ENV_R_PREFIX}"
```

## What each env is for

`pertussis-bio-tools`

- `Python 3.10`
- 常规生信工具和主分析依赖
- 包括 `blast`、`kraken2`、`mlst`、`checkm-genome`、`samtools`、`bcftools`、`bwa`、`freebayes`、`snippy`
- 也包括主分析脚本会直接用到的 `numpy`、`pandas`、`statsmodels`、`openpyxl`、`pycountry`、`snakemake`、`epydemix`

`pertussis-prn-global-bio`

- `Python 3.9`
- 主要承载 M4 系统发育和 M6 读段验证兼容层
- 当前配置除了 `gubbins`、`iqtree`、`raxml-ng`、`clonalframeml`、`pastml`、`ete3`、`ismapper` 之外，也补了 `panisa`、`emboss`、`samtools`、`bwa`、`bedtools`、`biopython`

`pertussis-prn-global-r`

- 独立 R 绘图与统计环境
- 包含 `tidyverse`、`ggplot2`、`ggrepel`、`patchwork`、`shiny`、`ape`、`ggtree`、`treeio`

## Quick checks

检查工具环境：

```bash
"${PROJECT_ENV_BIO_TOOLS_PREFIX}/bin/python" -c "import pysam, pandas, statsmodels, openpyxl, pycountry, epydemix"
"${PROJECT_ENV_BIO_TOOLS_PREFIX}/bin/snippy" --version >/dev/null
```

检查 3.9 兼容环境：

```bash
"${PROJECT_ENV_PHYLO_PREFIX}/bin/python" -c "import ete3, Bio"
"${PROJECT_ENV_PHYLO_PREFIX}/bin/iqtree2" -h >/dev/null
"${PROJECT_ENV_PHYLO_PREFIX}/bin/pastml" --help >/dev/null
```

检查 R 环境：

```bash
"${PROJECT_ENV_R_PREFIX}/bin/Rscript" -e "library(ggplot2); library(ggrepel); library(shiny); library(ggtree); library(treeio)"
```

## Note

Snakemake rules and the active shell wrappers now point directly at `config/env/` 下的这三份环境定义。
The older step-specific `envs/*.yml` specs were retired so dependency changes only need to be maintained in one place.

`checkm` 和 `gubbins` 当前都还会碰到 `pkg_resources` 兼容层，所以 `pertussis-bio-tools` 和 `pertussis-prn-global-bio` 都保留了 `setuptools<81` 这个兼容性约束。
