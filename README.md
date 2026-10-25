# PhyloMine

PhyloMine is a comprehensive comparative genomics and evolutionary biology workflow for ≤1000 MAGs or pure-culture genomes. It is specifically designed for deeply diverse datasets spanning classes, orders, and families, integrating a complete suite of quality control and dereplication, phylogenetic tree construction, and functional annotation to specialized functional mining (e.g., BGCs, CAZymes, defense systems) and HGT/DTL detection.

本流程基于 **Nextflow (DSL2)** 与 **Singularity** 容器化，覆盖从基因组质控、去冗余、物种树构建，到直系同源基因（OG）鉴定、功能注释、次级代谢基因簇挖掘、防御系统/CRISPR/HGT 检测，以及基因树与物种树 reconciliation 的完整分析链路。

---

## 目录

- [1. 分析流程](#1-分析流程)
- [2. 环境依赖](#2-环境依赖)
- [3. 数据库与配置文件](#3-数据库与配置文件)
- [4. 输入数据](#4-输入数据)
- [5. 快速开始](#5-快速开始)
- [6. 从中间步骤开始 / 断点续跑](#6-从中间步骤开始--断点续跑)
- [7. 参数说明](#7-参数说明)
- [8. 输出目录](#8-输出目录)
- [9. 模块与容器](#9-模块与容器)
- [10. 资源与 profile](#10-资源与-profile)
- [11. 常见问题与注意事项](#11-常见问题与注意事项)
- [12. 工具引用](#12-工具引用)

---

## 1. 分析流程

### 1.1 总体流程图

```
输入基因组目录 (--input)
        │
        ▼
┌─────────────────────────────────────────────────────────────────────┐
│ stage_1  质控、去冗余、物种树                                       │
│   CheckM2 → 质量过滤 → skDER 去冗余 → GTDB-Tk → GTDB 参考代表选择   │
│   → 物种树构建与 rooting                                            │
│   emit: rep_genomes, gtdbtk_out, species_root_tree,                 │
│         gtdb_rep_genomes_list, checkm2_report                       │
└─────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────┐
│ stage_2  基因预测与直系同源基因                                     │
│   Pyrodigal → OrthoFinder → OGFilter → MMseqs2 聚类                 │
│   emit: assets_ch, clean_OG, filtered_og_ch                         │
└─────────────────────────────────────────────────────────────────────┘
        │
        ├────────────────┬────────────────┬──────────────────────────┐
        ▼                ▼                ▼                          ▼
┌───────────────┐ ┌───────────────┐ ┌───────────────┐ ┌──────────────────────────┐
│ stage_3       │ │ stage_4       │ │ stage_5       │ │ stage_6                  │
│ OG 功能注释   │ │ BGC 挖掘      │ │ 防御系统/HGT  │ │ 基因树 + ALERAX          │
│ BAKTA         │ │ CGCFINDER     │ │ DefenseFinder │ │ IQ-TREE (per OG)         │
│ eggNOG        │ │ antiSMASH     │ │ CRISPRCasTyper│ │ ALERAX (DTL/transfer)    │
│ dbCAN         │ │ BiG-SCAPE     │ │ HGTector2     │ │                          │
│ RGI / VF      │ │               │ │               │ │                          │
│ MEROPS        │ │               │ │               │ │                          │
└───────────────┘ └───────────────┘ └───────────────┘ └──────────────────────────┘
```

### 1.2 各阶段说明

| Stage | 子流程 | 主要工作 |
|---|---|---|
| `stage_1` | `subworkflows/stage_1.nf` | CheckM2 质控，按完整度/污染度过滤，skDER 去冗余，GTDB-Tk 分类，按指定分类层级（`--level`）选择 GTDB 参考代表基因组，与用户基因组合并构建物种树并 rooting |
| `stage_2` | `subworkflows/stage_2.nf` | Pyrodigal 基因预测（meta 模式，输出 GFF、基因序列、蛋白序列）；OrthoFinder 鉴定直系同源组；`OGFilter.py` 保留单拷贝且覆盖 ≥8 个物种的 OG；MMseqs2 以 70% identity 聚类得到非冗余 OG 代表序列 |
| `stage_3` | `subworkflows/stage_3.nf` | 对 `clean_OG` 代表序列进行功能注释：BAKTA、eggNOG-mapper、dbCAN、RGI（CARD）、DIAMOND 比对 VFDB 与 MEROPS |
| `stage_4` | `subworkflows/stage_4.nf` | 对每个基因组进行次级代谢基因簇挖掘：dbCAN `cgc_finder`、antiSMASH，以及基于 antiSMASH 结果的 BiG-SCAPE 聚类 |
| `stage_5` | `subworkflows/stage_5.nf` | 对每个基因组进行防御系统（DefenseFinder）、CRISPR-Cas（CRISPRCasTyper）和 HGT（HGTector2）检测 |
| `stage_6` | `subworkflows/stage_6.nf` | 对每个 `filtered_OG` 用 MAFFT + IQ-TREE 建基因树（1000 次 UFBoot），再将基因树、物种树、物种缺失比例送入 ALERAX 进行 reconciliation / DTL 分析 |

> **注意**：`stage_3`、`stage_4`、`stage_5` 是 `stage_2` 的并列下游，三者互不依赖。`stage_4` 和 `stage_5` 需要的是 `stage_2` 的 `assets_ch`（每个基因组的 genome/faa/gff/genes），**不是** `total_OG_clean_70_rep_seq.fasta`。

---

## 2. 环境依赖

### 2.1 必需软件

- **Nextflow**（DSL2）：开发与测试版本为 `26.04.6`，建议使用相同或更新版本。
- **Singularity / Apptainer**：流程中的分析工具均通过容器运行。
  - 安装方法参考：<https://github.com/sylabs/singularity>
- **Slurm**（可选）：默认仍可通过 `-profile local` 在本机运行；集群运行使用 `-profile slurm`。
- **seqkit**：`GTDB_SELECT_REP` 直接在宿主机执行（该 process 未配置容器），需要 `seqkit` 在 `PATH` 中。

### 2.2 硬件建议

流程内存需求较高，单个任务峰值例如：

- GTDB-Tk：14 CPU / 120 GB
- ALERAX：12 CPU / 80 GB
- CheckM2：8 CPU / 40 GB
- skDER：8 CPU / 60 GB
- OrthoFinder：12 CPU / 60 GB
- Bakta / eggNOG：8 CPU / 60 GB

默认 `maxForks = 40`，实际部署时建议结合集群资源调整。

### 2.3 容器

`nextflow.config` 中已启用 Singularity：

```groovy
singularity {
    cacheDir = '/share/data01/project/xujinxin/project/20260910_nextflow_v2/container'
    enabled = true
    autoMounts = true
}
```

请把 `cacheDir` 改为本机/集群上可写、空间充足的目录。所有容器来自 `quay.io`；离线集群需要提前 `singularity pull` 到该缓存目录。

---

## 3. 数据库与配置文件

所有数据库路径在 `nextflow.config` 的 `params {}` 中配置。仓库中的默认值是开发环境的示例路径，**部署时必须改成自己的路径**。

| 参数 | 用途 | 挂载点 |
|---|---|---|
| `checkm2_db` | CheckM2 的 `.dmnd` 数据库文件 | `/opt/checkm2_db/uniref100.KO.1.dmnd` |
| `gtdbtk_db` | GTDB-Tk 参考数据库目录（示例为 GTDB R232） | `/opt/gtdbtk_db` |
| `gtdb_auxillary_metadata` | GTDB metadata TSV（用于选择参考代表） | 宿主机直接读取 |
| `gtdbtk_genome` | GTDB `genomic_fasta_reps` 目录 | 宿主机直接读取 |
| `bakta_db` | Bakta 数据库 | `/opt/bakta_db` |
| `eggNOG_db` | eggNOG-mapper 数据目录 | `/opt/eggNOG_db` |
| `dbcan_db` | dbCAN 数据库 | `/opt/dbcan_db` |
| `card_db` | CARD 数据库 | `/opt/card_db` |
| `vf_db` | VFDB DIAMOND 数据库文件 | `/opt/vf_db` |
| `antismash_db` | antiSMASH 数据库 | `/opt/antismash_db` |
| `pfam_db` | Pfam-A.hmm 所在目录（BiG-SCAPE 使用） | `/opt/pfam_db` |
| `macsyfinder_models` | DefenseFinder / MacSyFinder 模型目录 | `/opt/macsyfinder` |
| `cctyper_db` | CRISPRCasTyper 数据目录 | `/opt/cctyper_db` |
| `hgtector2_db` | HGTector2 数据库（含 `db/` 与 `taxdump/`） | `/opt/hgtector2_db` |
| `merops_db` | MEROPS DIAMOND 数据库文件 | `/opt/merops_db` |

各模块的容器与资源定义在 `modules/*.config` 中，由 `nextflow.config` 末尾的 `includeConfig` 统一加载。

---

## 4. 输入数据

- 一个基因组目录，包含待分析的基因组 FASTA 文件（`--input` 传入目录路径）。
- 默认扩展名为 `fa`（由 `genome_extension` 控制）。
- 每个文件一个基因组；**文件名（不含扩展名）作为样本 ID**，请勿包含空格。
- 输入应为**核苷酸基因组**，不是蛋白序列。
- 建议单次运行不超过 **1000** 个基因组。
- 所有基因组应来自**同一个域**（细菌或古菌）。`find_common_taxon.sh` 检测到 GTDB-Tk 同时产生 `bac120` 和 `ar53` 结果时会主动报错退出。

示例：

```bash
nextflow run main.nf \
    --input /path/to/genomes \
    --outdir /path/to/results \
    -profile slurm
```

---

## 5. 快速开始

### 5.1 完整流程（local）

```bash
nextflow run main.nf \
    --input /path/to/genomes \
    --outdir /path/to/results \
    -profile local
```

### 5.2 完整流程（Slurm）

```bash
nextflow run main.nf \
    --input /path/to/genomes \
    --outdir /path/to/results \
    -profile slurm
```

`slurm` profile 默认队列为 `gn`，请按集群实际情况修改 `nextflow.config` 中的 `queue`。

### 5.3 断点续跑

成功完成的任务会被缓存；使用相同的工作目录并加 `-resume`：

```bash
nextflow run main.nf \
    --input /path/to/genomes \
    --outdir /path/to/results \
    -profile slurm -resume
```

---

## 6. 从中间步骤开始 / 断点续跑

除了 `-resume`，流程内置了 `--run_stage` 参数，可以只执行某个阶段。当前支持 `stage_3` 和 `stage_4`。

### 6.1 只运行 stage_3

适用场景：已经有 `total_OG_clean_70_rep_seq.fasta`，只需要做功能注释。

```bash
nextflow run main.nf \
    --run_stage stage_3 \
    --clean_og /path/to/old_outdir/gencomp/total_OG_clean_70_rep_seq.fasta \
    --outdir /path/to/new_outdir \
    -profile slurm
```

前置条件：

- `clean_og` 指向 `stage_2` 输出的 `total_OG_clean_70_rep_seq.fasta`。

### 6.2 只运行 stage_4

适用场景：已经完成 `stage_2`，只需要做 BGC 挖掘。

```bash
nextflow run main.nf \
    --run_stage stage_4 \
    --assets_dir /path/to/old_outdir/gencomp \
    --outdir /path/to/new_outdir \
    -profile slurm
```

`assets_dir` 必须是上一次运行的 `gencomp` 目录，且包含：

```
gencomp/
├── rep_genomes/*.fa
└── pyrodigal/
    ├── proteins/*.faa
    ├── gff/*.gff
    └── genes/*.genes.fa
```

> 注意：`total_OG_clean_70_rep_seq.fasta` 是 `stage_3` 的输入，不是 `stage_4` 的输入。`stage_4` 需要的是每个基因组的 `genome/faa/gff/genes` 资产。

### 6.3 stage_5 / stage_6 单步运行

当前版本没有为 `stage_5`、`stage_6` 提供 `--run_stage` 入口。如需单独运行：

- 使用完整流程 + `-resume`，让已完成的阶段走缓存；
- 或者参考 `stage_4_only` 的写法，在 `main.nf` 中新增 `stage_5_only` / `stage_6_only` 分支。

`stage_6`（IQ-TREE + ALERAX）除 `filtered_og_ch` 外，还需要 `stage_1` 的：

- `species_root_tree`（`gtdbtk.infer/gtdbtk.rooted.tree`）
- `gtdb_rep_genomes_list`（`gtdb_rep_genomes/representative_genome.tsv`）
- `checkm2_report`（`checkm2/quality_report.tsv`）

### 6.4 关于 `publishDir` 和软链接

所有结果通过 `publishDir "${params.outdir}/gencomp/", mode: 'symlink'` 发布。`gencomp/` 中的文件可能是指向 `work/` 的软链接：

- 保留 `work/` 目录，否则旧结果的软链接会失效；
- 需要长期归档时，请复制实体文件（如 `cp -L`）。

---

## 7. 参数说明

### 7.1 核心参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--input` | `null` | 输入基因组目录 |
| `--outdir` | `null` | 输出根目录；结果写入 `${outdir}/gencomp` |
| `--genome_extension` | `fa` | 基因组文件扩展名 |
| `--min_completeness` | `70` | CheckM2 完整度阈值（%） |
| `--max_contamination` | `5` | CheckM2 污染度阈值（%） |
| `--level` | `order` | GTDB 参考代表选择的分类层级：`phylum` / `class` / `order` / `family` / `genus` / `species` |
| `--run_stage` | `null` | 只运行单个阶段：`stage_3` 或 `stage_4` |
| `--clean_og` | `null` | `stage_3` 单步运行时使用的 `total_OG_clean_70_rep_seq.fasta` |
| `--assets_dir` | `null` | `stage_4` 单步运行时使用的上次 `gencomp` 目录 |

### 7.2 数据库参数

见 [3. 数据库与配置文件](#3-数据库与配置文件)。

### 7.3 流程行为参数

以下参数在 `nextflow.config` 中：

| 配置 | 默认值 | 说明 |
|---|---|---|
| `process.executor` | `local` | 执行器，`slurm` profile 会覆盖 |
| `process.cpus` | `4` | 默认 CPU |
| `process.memory` | `20.GB` | 默认内存 |
| `process.time` | `24.h` | 默认运行时间 |
| `process.maxForks` | `40` | 最大并行任务数 |
| `process.errorStrategy` | `finish` | 出错后等待已启动任务结束再退出 |

各 process 在 `modules/*.config` 中通过 `withLabel` 覆盖 CPU/内存/容器。

---

## 8. 输出目录

所有结果统一发布到 `${params.outdir}/gencomp/`：

```
<outdir>/
└── gencomp/
    ├── checkm2/
    │   ├── quality_report.tsv
    │   └── ...
    ├── high_quality_genomes/
    ├── skder/
    │   └── Dereplicated_Representative_Genomes/
    ├── gtdbtk/
    ├── rep_genomes/                        # 下游使用的基因组（用户 + GTDB 参考）
    ├── gtdb_rep_genomes/
    │   └── representative_genome.tsv
    ├── gtdbtk.infer/
    │   └── gtdbtk.rooted.tree              # species_root_tree（树标签已去除 RS_/GB_ 前缀）
    ├── pyrodigal/
    │   ├── gff/<id>.gff
    │   ├── genes/<id>.genes.fa
    │   └── proteins/<id>.faa
    ├── orthofinder/
    │   ├── OGfilter/OG*.fa
    │   └── ...
    ├── total_OG.clean.fa
    ├── total_OG_clean_70_rep_seq.fasta     # clean_OG
    ├── bakta/
    ├── eggNOG/
    ├── dbcan/
    ├── card/
    ├── vf/
    ├── merops/
    ├── cgc/<genome>/
    ├── antismash/<genome>/
    ├── bigscape/
    ├── defensefinder/<genome>/
    ├── crisprcastyper/<genome>/
    ├── hgtector2/<genome>/
    ├── iqtree_og/<OG>/
    ├── family.txt                          # ALERAX 输入
    ├── missing.txt                         # ALERAX 输入的物种缺失比例
    └── output_3_per_family*                # ALERAX reconciliation 结果
```

---

## 9. 模块与容器

| Stage | Process | 工具 / 版本 | 容器 |
|---|---|---|---|
| 1 | `CHECKM2` | CheckM2 1.1.0 | `quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1` |
| 1 | `FILTER_GENOMES` | awk / cp | 宿主机 |
| 1 | `SKDER` | skDER 1.3.7 | `quay.io/biocontainers/skder:1.3.7--py310hd766df8_0` |
| 1 | `GTDBTK` | GTDB-Tk 2.7.2 | `quay.io/biocontainers/gtdbtk:2.7.2--pyhdfd78af_0` |
| 1 | `GTDB_SELECT_REP` | find_common_taxon / fetch_rep_genome / seqkit | 宿主机 |
| 1 | `ROOT_TREE` | GTDB-Tk 2.7.2 | `quay.io/biocontainers/gtdbtk:2.7.2--pyhdfd78af_0` |
| 2 | `PYRODIGAL` | Pyrodigal 3.7.1 | `quay.io/biocontainers/pyrodigal:3.7.1--py310h7c10099_1` |
| 2 | `ORTHOFINDER` | OrthoFinder 2.5.5 | `quay.io/biocontainers/orthofinder:2.5.5--hdfd78af_2` |
| 2 | `CLUSTER_ORTHOLOGS` | MMseqs2 | `quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0` |
| 3 | `BAKTA` | Bakta 1.12.1 | `quay.io/biocontainers/bakta:1.12.1--pyhdfd78af_0` |
| 3 | `EGGNOG` | eggNOG-mapper 3.0.0-beta6 | `quay.io/biocontainers/eggnog-mapper:3.0.0-beta6` |
| 3 | `CAZYME` | dbCAN 5.2.9 | `quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0` |
| 3 | `RGI` | RGI 6.0.8 | `quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0` |
| 3 | `VF` | DIAMOND 2.2.6 | `quay.io/biocontainers/diamond:2.2.6--he361c42_0` |
| 3 | `MEROPS` | DIAMOND 2.2.6 | `quay.io/biocontainers/diamond:2.2.6--he361c42_0` |
| 4 | `CGCFINDER` | dbCAN 5.2.9 | `quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0` |
| 4 | `ANTISMASH` | antiSMASH 8.0.4 | `quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_0` |
| 4 | `BIGSCAPE` | BiG-SCAPE 2.0.3 | `quay.io/biocontainers/bigscape:2.0.3--pyhdfd78af_0` |
| 5 | `DEFENSEFINDER` | DefenseFinder 3.0.0 | `quay.io/biocontainers/defense-finder:3.0.0--pyhdfd78af_0` |
| 5 | `CCTYPER` | CRISPRCasTyper 1.9.0 | `quay.io/mahoon/cctyper:1.9.0` |
| 5 | `HGTECTOR2` | HGTector2 0.2.1 | `quay.io/mahoon/hgtector:0.2.1` |
| 6 | `IQTREE_OG` | MAFFT + IQ-TREE（OrthoFinder 容器） | `quay.io/biocontainers/orthofinder:2.5.5--hdfd78af_2` |
| 6 | `ALERAX` | ALERAX 1.4.0 | `quay.io/biocontainers/alerax:1.4.0--heeb98d5_0` |

> `CLUSTER_ORTHOLOGS` 在 `modules/orthofinder.nf` 中定义，但容器配置位于 `modules/orthofinder.config`。

---

## 10. 资源与 profile

### 10.1 内置 profile

| Profile | 说明 |
|---|---|
| `local` | 本机执行 |
| `slurm` | Slurm 执行；默认队列 `gn`，`queueSize=100`，`submitRateLimit=10/1sec` |

### 10.2 资源覆盖

`nextflow.config` 的默认值为 4 CPU / 20 GB / 24 h，各 process 的 `withLabel` 配置会覆盖。若集群资源与默认值不同，建议：

- 修改 `modules/*.config` 中的 `cpus` / `memory`；
- 或使用自定义 config 覆盖，例如：

```groovy
process {
    withLabel: 'gtdbtk' { memory = '100 GB' }
    withLabel: 'alerax' { cpus = 8; memory = '60 GB' }
}
```

运行：

```bash
nextflow run main.nf -c my_resources.config ...
```

### 10.3 ALERAX 的 MPI 参数

ALERAX 通过 `mpiexec -np ${task.cpus}` 启动。如果 OpenMPI 报：

```
There are not enough slots available in the system ...
```

说明容器内的 OpenMPI 没有识别到 Slurm 分配，按物理核数计算 slot。流程中已使用：

```bash
mpiexec --use-hwthread-cpus --oversubscribe -np ${task.cpus}
```

若仍有问题，可按集群情况调整 `cpus` 或 MPI 参数。

---

## 11. 常见问题与注意事项

### 11.1 部署前必须修改的配置

1. `nextflow.config` 中所有数据库路径；
2. `singularity.cacheDir`；
3. Slurm profile 的 `queue`（默认 `gn`）；
4. 视集群资源调整各模块的 `cpus` / `memory`。

### 11.2 基因组 ID 与 `RS_` / `GB_` 前缀

GTDB 参考基因组在 metadata / GTDB-Tk 结果中带有 `RS_`（representative species）或 `GB_`（GenBank representative）前缀，但真实的 NCBI accession 是 `GCF_...` / `GCA_...`。

流程的处理方式：

- `fetch_rep_genome.sh` 在复制 GTDB 参考基因组时去掉前缀，`rep_genomes/*.fa` 使用 `GCF_...` / `GCA_...`；
- `ROOT_TREE` 在生成 `gtdbtk.rooted.tree` 后会把树标签中的 `RS_` / `GB_` 去掉；
- `ALERAX` 的 `mapping` 和 `missing.txt` 同样使用去前缀后的 ID。

因此三处物种 ID 保持一致。升级或替换 GTDB 数据后，建议先检查去前缀是否产生重复：

```bash
cut -f1 gencomp/gtdb_rep_genomes/representative_genome.tsv \
  | sed -E 's/^(RS|GB)_//' \
  | sort | uniq -d
```

无输出表示没有重复。

### 11.3 `-stub` 模式的限制

仓库中部分 process 的 `stub` 块不完整，例如：

- `GTDB_SELECT_REP` 的 stub 只创建 `gtdbtk`，但声明的输出是 `rep_genomes/*` 和 `gtdb_rep_genomes/representative_genome.tsv`；
- `RGI` 的 stub 创建 `rgi`，但输出声明是 `card`；
- `IQTREE_OG` 的 stub 只创建 `iqtree_og`，没有 `iqtree_og/${id}/${id}.ufboot`。

因此 `-stub-run` 可能在中途报 `Missing output file(s)`。这只是测试桩的问题，不影响真实数据运行；不要用 stub 结果判断流程是否正确。

### 11.4 常见运行错误

| 现象 | 原因 | 处理 |
|---|---|---|
| `Invalid process definition -- check for missing or out-of-order section labels` | process 内有空的 `output:` 段，或 section 顺序错误 | 确保顺序为 `input:` → `output:` → `script:` → `stub:`，且每段非空 |
| `No such property: xxx for class: ...ChannelOut` | `main.nf` 引用了子流程 `emit:` 中不存在的输出 | 检查对应 `subworkflows/stage_x.nf` 的 `emit:` 块 |
| OpenMPI `not enough slots` | 容器内 OpenMPI 未识别 Slurm 分配 | 使用 `--use-hwthread-cpus --oversubscribe`（已内置） |
| `publishDir` 结果链接失效 | `work/` 被清理 | 保留 `work/`，或用 `cp -L` 归档实体文件 |
| GTDB-Tk 同时出现 `bac120` 与 `ar53` | 输入混合了细菌和古菌 | 分开运行，单次只处理一个域 |

### 11.5 中间结果复用

- `stage_3` 只需要 `total_OG_clean_70_rep_seq.fasta`；
- `stage_4` 只需要 `rep_genomes/` 和 `pyrodigal/{proteins,gff,genes}/`；
- `stage_6` 需要 `stage_1` 的物种树、CheckM2 报告、GTDB 代表列表，以及 `stage_2` 的 `filtered_OG` 与 IQ-TREE 基因树；
- 所有阶段都建议配合 `-resume` 使用，避免重复计算。

---

## 12. 工具引用

使用本流程时，请引用其中使用的分析工具。各工具的官方项目与引用信息见：

| 工具 | 项目地址 |
|---|---|
| Nextflow | https://www.nextflow.io |
| CheckM2 | https://github.com/chklovski/CheckM2 |
| GTDB-Tk | https://github.com/Ecogenomics/GTDBTk |
| skDER | https://github.com/raufs/skDER |
| Pyrodigal | https://github.com/althonos/pyrodigal |
| OrthoFinder | https://github.com/davidemms/OrthoFinder |
| MMseqs2 | https://github.com/soedinglab/MMseqs2 |
| Bakta | https://github.com/oschwengers/bakta |
| eggNOG-mapper | https://github.com/eggnogdb/eggnog-mapper |
| dbCAN | https://github.com/linnabrown/run_dbcan |
| RGI / CARD | https://github.com/arpcard/rgi |
| antiSMASH | https://github.com/antismash/antismash |
| BiG-SCAPE | https://github.com/medema-group/BiG-SCAPE |
| DefenseFinder | https://github.com/mdmparis/defense-finder |
| CRISPRCasTyper | https://github.com/Russel88/CRISPRCasTyper |
| HGTector | https://github.com/qiyunzhu/HGTector |
| IQ-TREE | https://github.com/iqtree/iqtree2 |
| MAFFT | https://github.com/GSLBiotech/mafft |
| ALERAX / ALE | https://github.com/maxemil/ALE |

请以各工具官方文档给出的最新文献为准。
