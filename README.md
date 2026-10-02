# CysRank: Structure-Aware Cysteine Site Prioritization

当前正式模型为 **sst-ranking-20260924-v1**：保留 22 个分类监督位点，新增一条 SST 实验排序监督。SST 已属于训练数据，99/100 构象排序符合仅为训练回代；见 [当前训练报告](tmp/docs/SST_TRAINING.md)。下文早期 MD/PROPKA 对照记录属于旧版模型。

面向多肽环化研究，本项目将蛋白质预训练模型的结构表征与位点生化特征相结合，建立从实验数据组织、模型训练到候选排序和结果追溯的计算流程，为后续实验提供可复核的半胱氨酸位点优先级参考。

系统支持对已有半胱氨酸（Cys）位点评分，也支持在输入窗口提供陌生多肽后，自动枚举新增一个 Cys 的候选序列、调用本地 ESMFold 预测结构，并输出 Top3 新序列及外部能垒评估交接包。通过过滤的候选不足三条时，按内部评分补足，并提示成功可能性较低及原因。另提供外部模拟能垒反馈与独立回归接口。

## 陌生序列插入设计

**前端使用说明已保存在 [自动结构预测与多肽输入窗口](tmp/docs/AUTOMATIC_DESIGN.md)。** 该说明包含环境安装、启动参数、输入限制、排名规则和外部评估包内容；按当前赛事目录要求保留在 `tmp/docs/`，不在项目根目录创建 `docs/`。

网页操作步骤：

1. 在已安装 `requirements.txt` 的项目 Python 环境中执行 `python -m pip install -r requirements-folding.txt`，然后运行 `python app.py`。保持服务终端运行，在浏览器打开 [本地操作界面](http://127.0.0.1:8765)。首次使用需要下载 ESMFold 权重，缓存完整后可运行 `python app.py --local-files-only`。
2. 在“目标多肽序列”中粘贴一条序列，例如 `AGKSTV`。只输入 20 种标准氨基酸单字母代码，不含 FASTA 标题、编号或端基修饰标记；空格和换行自动忽略，小写自动转大写。默认接受 2–100 个残基。
3. 点击“开始预测”，等待结构预测和评分完成。可展开“查看运行详情”检查进度；同一服务一次只处理一个任务，CPU 预测可能较慢。连接中断时点击“重新连接任务”，服务仍在运行时也可刷新页面恢复当前任务。
4. 在“推荐序列 / TOP 3”中查看候选。高亮的 `C` 是新增残基，候选长度比输入多 1；“新增 Cys”编号从候选序列第 1 位开始计数。优先展示通过过滤的候选，不足时按分数补足，并显示原因；不同序列不足三条时显示实际数量。模型分数用于排序，不是实验成功概率，本设计流程不以 0.5 为入选硬门槛。
5. 点击“下载结构”保存候选 PDB，点击“查看报告”下载 Markdown 报告，或点击“下载外部评估包”保存包含 Top3 CSV、结构、评分和交接文件的 ZIP。网页不执行 GFN2-xTB 能垒计算。
6. 默认结果保存在 `results/web_jobs/<任务编号>/design/`，日志在同一任务目录的 `progress.log`。服务重启后文件仍在，但网页不能恢复旧任务；请从磁盘读取结果。若显示“预测未完成”，先检查日志中的依赖、模型下载或内存错误，处理后重新提交。

**当前目录迁移的已知影响：** `src/folding.py` 导出评估包时仍读取根目录下的 `docs/BARRIER_INTERFACE.md`，文件迁移后会导致导出阶段失败；`tools/package_submission.py` 也仍要求根目录存在 `docs/`。运行完整预测或打包前，需要将这些代码中的旧路径与当前目录布局同步。本次更新为文档入口整理，尚未修改这两处代码，也未重新执行完整预测。

```powershell
python design.py run --sequence AGKSTV --record-id demo_only --output results/auto_design_demo
```

以下手工提供结构的流程仍可使用：

```powershell
python design.py prepare --sequence AGKSTV --record-id demo_only --output results/design_demo_requests
# 填写生成的结构清单，提供每条插入后候选的真实配套 PDB
python design.py score --candidates results/design_demo_requests/candidates.jsonl --structures results/design_demo_requests/structures.json --output results/design_demo_scored
```

上述序列仅用于流程示例。新 Cys 为新增残基，输出长度增加 1；原序列 PDB 不能直接代替插入后结构。评分后导出 Top3 CSV、结构证据、PDB 副本和待填写的能垒请求。完整使用方法及结构来源要求见 [插入设计说明](tmp/docs/DESIGN_WORKFLOW.md)。下文 `predict.py` 和 `run.py` 仍是已有位点评分与训练复现入口。

## 项目特色

- **融合序列、结构与生化信息。** 同时利用 ProteinMPNN 结构嵌入、局部序列上下文，以及二硫键状态、溶剂可及表面积、前位脯氨酸和邻近原子数，为位点评分提供互补信息。
- **轻量化任务适配。** 冻结预训练编码器，仅训练包含 364 个参数的分类头，在复用结构表征的同时控制可训练参数规模，便于在有限标注数据下开展方法探索。
- **预测与结构约束协同筛选。** 将模型分数与显式过滤规则结合，输出候选排序，同时保留每个位点的评分、过滤状态和拒绝原因，便于人工复核。
- **区分实验结局与未验证位点。** 使用独立监督掩码，避免将未标注 Cys 直接视为阴性；反馈接口区分湿实验、模拟与留出数据，支持有来源记录的迭代训练。
- **完整的可复现交付。** 提供训练入口、配置、最终权重、标准化结果文件和验证记录；所有运行相对路径均以仓库根目录为基准，便于迁移到评审环境。

## 技术方案

```text
多肽序列 + Cα 坐标 + 逐残基生化特征
                    │
       ┌────────────┼────────────┐
       │            │            │
 ProteinMPNN    局部序列窗口    四项生化特征
 结构嵌入       11 个残基       逐位点输入
 128 维          231 维          4 维
       └────────────┼────────────┘
                    │
             363 维特征融合
                    │
            线性分类头 + Sigmoid
                    │
         固定阈值 + 二硫键 / 前位 Pro 过滤
                    │
         Top 3 候选清单 + 全位点评分记录
```

### 结构表征与任务学习

项目采用 CA-only ProteinMPNN 作为冻结的结构编码器，并在其上训练面向 Cys 位点的二分类头。结构嵌入描述残基周围的空间环境，11 残基窗口编码局部序列上下文，四项生化特征补充位点的可及性与局部环境信息。

训练使用实验来源记录中明确标注的位点，标签和样本来源不作为预测特征。模型采用 AdamW 优化器，固定随机种子、训练轮数和判定阈值，并提供按母蛋白分组的留出诊断和特征消融结果。

### 候选筛选与结果追溯

推理阶段先对指定 Cys 位点评分，再依据二硫键状态和前位 Pro 标记进行过滤。通过过滤且分数不低于 0.5 的位点按分数降序排列，默认输出全局 Top 3；同分时按记录编号及位点位置排序。前位 Pro 规则作为当前采用的保守筛选假设，SASA 与邻近原子数参与模型输入，不设置额外硬阈值。

每次运行同时保存候选清单、全部位点评分和运行元数据。候选可追溯至输入记录、模型版本及权重哈希，便于比较不同运行结果。

### 实验反馈接口

系统支持对后续实验或模拟观察进行格式校验与归档。只有经人工审核、归属于训练集的湿实验反馈可用于迭代训练；模拟、留出和未审核记录独立保留。重复观察、重复序列位点和留出母蛋白冲突会被检查，为后续实验反馈接入提供基础。

## 快速开始

### 1. 配置环境

已验证环境如下；完整依赖版本见 [requirements.txt](requirements.txt)。

| 项目 | 配置 |
| --- | --- |
| 操作系统 | Windows 11 x64，build 26200 |
| Python | 3.13.5 |
| PyTorch | 2.11.0+cu128 |
| NumPy | 2.5.2 |
| GPU | NVIDIA GeForce RTX 4050 Laptop GPU，6 GB |
| NVIDIA 驱动 | 592.82 |
| CUDA 运行时 | PyTorch CUDA 12.8 |

在仓库根目录创建独立环境：

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements.txt
```

以下命令中的 `python` 指已安装项目依赖的解释器，例如 Windows 下的 `.venv/Scripts/python.exe`。推理使用 CPU，训练自动选择可用的 CUDA 或 CPU，分类头在 CPU 上优化。项目运行不依赖在线 API，也无需调用 PyMOL 或 ChimeraX；输入结构特征需预先准备。

依赖文件锁定上述已验证环境。其他平台使用对应的 PyTorch 分发包后需重新核验运行结果；Linux 尚未实机验证。

### 2. 生成候选清单

```bash
python predict.py --input data/examples/input.jsonl --output results/review.csv --track "实际参赛赛道"
```

运行后生成：

| 文件 | 内容 |
| --- | --- |
| `results/review.csv` | Top 3 候选位点、预测分数、输入结构特征及版本信息 |
| `results/review.sites.jsonl` | 全部已评分位点、分类结果和过滤原因 |
| `results/review.run.json` | 输入与模型哈希、软件版本、筛选阈值及候选数量 |

项目已提供 [示例候选清单](results/results.csv)，可直接查看输出格式。示例输入来自训练数据，用于演示运行流程。

### 3. 一键复现完整流程

```bash
python run.py --output logs/demo --track "实际参赛赛道"
```

该命令完成 **模型训练 → 参数优化 → 位点评分 → 候选筛选**，保存训练权重、损失日志、分组评估、候选清单与运行耗时。默认运行 50 轮，使用项目随附的小规模数据；`--epochs 1` 可用于快速检查环境。

**路径约定：所有配置文件及命令行中的相对路径均基于仓库根目录解析，与终端当前工作目录无关。** 从其他位置运行时，只需正确指定入口脚本的位置。

为保护已有结果，输出文件或训练输出目录须为新路径。默认预测位置为 `results/results.csv`，该文件已随项目提供，复跑时请使用上例的新输出路径。`--track` 应填写实际参赛赛道；未填写时标记为演示。

## Jupyter Notebook

[打开完整流程 Notebook](notebooks/01_reproducible_workflow.ipynb)。包含中文步骤说明、数据核对、50 轮训练与参数优化、母蛋白分组诊断、随附最终模型与新训练模型的 Top 3 对照，以及原始评分一致性核验。已保存实际执行输出。

使用项目统一依赖文件配置环境并启动：

```bash
python -m pip install -r requirements.txt
python -m jupyterlab
```

选择与项目依赖相同的 Python 内核，打开 Notebook 后执行 **Restart Kernel and Run All Cells**。默认 50 轮，参数单元格可改为 1 轮冒烟检查。可从仓库根目录或 notebooks 目录启动；内核工作目录在仓库之外时，先设置 `CYSRANK_ROOT` 指向移植后的仓库。

所有数据、模型和输出路径仍基于仓库根目录。每次完整执行都会写入新的 `logs/notebook/<运行编号>/`，不会替换正式权重或覆盖已有候选清单。Notebook 中的候选属于训练回代示例，不能当作独立验证。旧 notebook 与旧输出已归档至 tmp/archive/before_part2_retrain_20260924；当前 notebook 保存新一轮执行输出。

## 训练与迭代

### 独立训练

```bash
python train.py --config configs/train.json --output logs/reproduction
```

| 参数 | 默认设置 |
| --- | --- |
| 编码器 | 冻结的 CA-only ProteinMPNN |
| 分类头 | 363 维输入，364 个可训练参数 |
| 优化器 | AdamW |
| 学习率 / 权重衰减 | 0.01 / 0.01 |
| 训练轮数 | 50 |
| 随机种子 | 20260916 |
| 判定阈值 | 0.5 |

训练配置位于 [configs/train.json](configs/train.json)。最终模型使用固定训练轮次，不根据留出诊断结果调整阈值或选择权重。

### 独立筛选

```bash
python screen.py --input data/examples/input.jsonl --output results/screen.csv --track "实际参赛赛道" --top-k 3
```

`screen.py` 与 `predict.py` 使用一致的筛选逻辑，支持通过 `--top-k` 控制候选数量。候选不足时输出实际符合条件的位点。

### 接入反馈

```bash
python feedback.py --input data/feedback/incoming.jsonl --output data/feedback/reviewed.jsonl
python train.py --feedback data/feedback/reviewed.jsonl --output logs/feedback_iteration
```

反馈文件由使用者提供，字段与审核要求见 [反馈接口说明](tmp/docs/FEEDBACK.md)。

### 外部能垒模拟反馈

已提供独立的模拟能垒训练接口：`barrier.py requests` 导出计算请求，计算人员交回结果后由 `feedback.py` 校验，`barrier.py train` 训练能垒回归器，`barrier.py predict` 为同一计算协议下的新位点评估能垒。支持 GFN2-xTB 或 DFT 结果，电子能垒与活化自由能分别建模。模拟不进入原有实验分类头，计算失败和未审核结果不训练。交接字段、命令和限制见 [能垒接口说明](tmp/docs/BARRIER_INTERFACE.md)。量化计算由外部执行，本项目尚无真实能垒反馈或已验证的能垒模型。

## 数据与结果格式

### 输入

输入采用 UTF-8 JSONL 格式，每行描述一条序列。设序列长度为 N，序列至少包含 3 个残基。

| 字段 | 形式 | 含义 |
| --- | --- | --- |
| `record_id` | 唯一非空字符串 | 输入记录编号 |
| `seq` | 长度 N 的字符串 | 氨基酸序列 |
| `xyz` | N × 3 数组 | Cα 坐标，单位 Å |
| `is_disulfide` | N 个 0/1 | 二硫键状态 |
| `sasa` | N 个非负数值 | 溶剂可及表面积，单位 Å² |
| `prev_is_pro` | N 个 0/1 | 前一残基是否为 Pro |
| `nearby_atoms` | N 个非负数值 | 5 Å 范围内非氢原子数 |
| `feature_mask` | N 个 0/1 | 标记结构特征已准备、需评分的位点 |
| `structure_file` | 可选相对路径 | 仓库内对应的结构文件 |

模型只评分 `feature_mask=1` 的 Cys 位点，推理不需要标签。SASA 为 0 可以是有效测量；缺失特征应补齐后再输入，不能用 0 代替。结构文件应与输入记录真实对应。

### 输出

候选 CSV 包含候选编号、所属赛道、序列、Cys 位置、关键预测分数、排序、模型版本、运行版本、权重哈希、结构文件引用、四项结构特征和备注。位点采用 1-based 编号，分数保留 8 位小数。

评分用于候选优先级比较。`predict.py` 评估已有 Cys 位点；`design.py run` 自动生成候选配套结构并完成评分与筛选，`design.py score` 支持手工结构输入。GFN2-xTB 计算及实验验证由外部完成。

## 实验结果与可复现性

本轮完成正例 SASA 特征更新后的数据整合，采用 19 条记录、22 个监督位点（8 正、14 负）训练轻量分类头。固定 50 轮训练后，BCE 从 **0.592711 降至 0.051014**，训练回代 **22/22** 正确；母蛋白分组诊断命中 **3/5 个可评估正位点**。

已完成的工作包括：

- **数据管理：** 统一监督标签语义，保留清洗记录、来源哈希和版本快照。
- **模型训练：** 融合预训练结构表征、序列上下文和生化特征，以 364 个可训练参数完成任务适配。
- **候选筛选：** 统一评分、规则过滤与 Top 3 导出，同时保存全部位点评分和筛选依据。
- **可复现验证：** 28 项自动测试通过；完整模型 CPU 重载通过，Notebook 从头执行后与正式训练指标和分组诊断一致。
- **拓展接口：** 支持插入候选枚举、外部配套结构评分及实验反馈数据接入。

评估范围：22/22 为训练拟合结果，3/5 为小样本正例分组诊断，无独立阴性测试；两者不代表新游离短肽的总体准确率。本轮同时更新特征与样本组成，不将指标差异单独归因于 SASA。适用范围见 Model Card。

```bash
python -m unittest discover -s tests -v
```

完整说明见 [Model Card](models/MODEL_CARD.md)、[训练报告](logs/retrain_20260924/TRAINING_REPORT.md) 和 [验证记录](logs/verification_20260924/verification.json)。

## 项目结构

```text
CysRank-GCLSC/
├── README.md                 # 项目介绍与使用指南
├── requirements.txt          # 运行依赖及版本
├── configs/                  # 训练配置
├── data/                     # 原始数据、处理后数据、示例及结构文件
├── src/                      # 特征提取、模型训练与推理逻辑
├── models/                   # 最终模型、预训练权重及 Model Card
├── train.py                  # 训练入口
├── predict.py                # 标准化结果生成入口
├── screen.py                 # 候选筛选入口
├── run.py                    # 完整流程入口
├── design.py                 # 单 Cys 插入设计及结构评分入口
├── barrier.py                # 能垒请求、训练和预测入口
├── feedback.py               # 反馈校验与归档入口
├── results/                  # 示例结果与位点评分
├── logs/                     # 训练、评估与验证记录
├── tests/                    # 自动化测试
├── tmp/docs/                 # 本地使用说明与参考资料，不纳入提交包
├── notebooks/               # 可执行的训练与筛选演示
├── tools/                   # 本地提交包生成工具
└── licenses/                 # 第三方许可证
```

## 还原态测试肽：Cys3 与 Cys14

已按提供方确认的 Ac/NH2 端基、还原巯基、pH 7.5、300 K，完成 10 起点短时显式水模拟。100 个生产快照中仅 7 个为 Cys3 更高，平均分 Cys3=0.8918、Cys14=0.9853；按起点平均，10 条均为 Cys14 更高。因此当前计算未稳定复现实验排序，不能称为独立有效性验证。

[查看案例与复现命令](tmp/docs/REDUCED_PEPTIDE_CASE.md) · [打开案例 Notebook](notebooks/02_reduced_peptide_case.ipynb) · [逐位点结果与结构关联](results/2mi1_reduced_capped_ph75_20260924/analysis/results.csv)

只复算随附结构：`python tools/analyze_2mi1_reduced_md.py --output results/2mi1_reanalysis`。重做模拟另安装 `requirements-md.txt`；详见案例指南中的 OpenCL 配置。该案例不插入新残基、不预测反应能垒，不改变正式模型。

## 本地整理与打包

专题说明现位于 `tmp/docs/`，README 中的相关链接已指向该位置；[历史整理说明](tmp/docs/PROJECT_LAYOUT.md) 中关于根目录 `docs/` 的描述反映迁移前布局。重复运行产物、历史版本、编辑器配置和过程文件统一归档到 tmp。环境按 requirements.txt 配置；旧路径对预测导出和打包的影响见上文。

另提供 [KTTKS 文献案例](tmp/docs/KTTKS_CASE_STUDY.md)及 results/kttks_design_requests/ 下的六条待评估序列，目前没有候选结构或排名。

```bash
python tools/package_submission.py --output tmp/agents/submission/CysRank-GCLSC.zip
```

打包工具使用明确的文件目录范围并生成逐文件 SHA-256 清单，排除 tmp、.git、本机环境、缓存及重复 Notebook 运行目录。只生成本地检查包，不上传；输出须为新路径。历史完整流程耗时约 12–37 秒，实际耗时随硬件和加载状态变化。

## 数据来源与第三方声明

项目基于 ProteinMPNN 的预训练结构编码能力，增加 Cys 位点监督分类头、生化特征融合、候选筛选及反馈管理流程。ProteinMPNN 的 MIT 许可保留于 [licenses/ProteinMPNN-MIT.txt](licenses/ProteinMPNN-MIT.txt)，本地模型版本与文件哈希见 [第三方说明](tmp/docs/THIRD_PARTY.md)。

数据来源、清洗、划分及许可信息见 [数据说明](data/README.md)。部分原始实验来源、数据授权和获取时间尚需数据提供方补充；正式提交范围及候选字段需与实际赛道要求核对，详见 [提交要求核对表](tmp/docs/SUBMISSION_CHECKLIST.md)。`tmp/`、`.git/`、本机虚拟环境与缓存不属于提交内容。


PROPKA 固定 λ=0.5 的探索性对照已完成：Cys3 更高的生产快照由 7/100 变为 6/100，未改善实验排序。包含端基类型适配及其局限，见 [完整报告](results/2mi1_propka_lambda05_20260924/REPORT.md)。


## 当前模型：SST 排序监督版

正式权重已更新为 sst-ranking-20260924-v1。原 22 个分类标签不变，新增一条 Cys3 > Cys14 实验排序监督；100 个构象共同占一条观察的权重。SST 训练回代 99/100 帧排序符合，不能作为独立验证。此前还原态/PROPKA 报告保留为旧模型对照。见 [重训说明](tmp/docs/SST_TRAINING.md)。
