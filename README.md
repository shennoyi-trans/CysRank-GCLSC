# CysRank｜面向多肽环化的结构感知半胱氨酸位点优选

**从一条多肽序列出发，生成 Cys 插入候选，结合三维结构与生化特征给出实验优先级，并保留可复核的筛选依据。**

CysRank 面向多肽环化研究中的候选选择问题：已有多个半胱氨酸（Cys）时，哪些位点值得优先研究？需要引入新的 Cys 时，应优先尝试哪些插入位置？项目将预训练蛋白质结构表征、轻量监督学习、显式结构检查和实验反馈接口整合为一套可运行、可追溯的辅助筛选流程，为后续实验与外部计算提供候选清单和结构材料。

项目支持两类任务：**对已有 Cys 位点评分排序**，以及**为新输入的多肽枚举单 Cys 插入方案并推荐 Top 3 序列**。使用者既可以通过网页输入序列，也可以通过命令行、JSONL 数据和 Jupyter Notebook 完成批量评分与训练复现。

当前正式模型为 **`sst-ranking-20260924-v1`**，在 22 个分类监督位点之外引入一条 SST 实验排序观察。模型分数表示候选优先级，不是实验成功概率；当前结果的验证范围见下文“成果与验证”。

**阅读导航：** [核心特色](#核心特色) · [成果与验证](#成果与验证) · [技术方案](#技术方案) · [快速开始与评审演示](#快速开始与评审演示) · [数据与结果格式](#数据与结果格式)

## 项目要解决什么问题

围绕“先做哪些实验”，CysRank 将候选生成、位点评分和结果交接连接起来，让研究者能够从候选序列进一步查看结构、筛选依据和数据来源。

| 研究与使用需求 | CysRank 的实现 | 对研究流程的价值 |
| --- | --- | --- |
| 同一多肽存在多个潜在插入位置 | 对长度为 N 的序列枚举 N+1 个单 Cys 插入操作，保留位置映射 | 系统覆盖候选位置，减少手工枚举与编号工作 |
| 位点选择需要考虑局部空间环境 | 融合结构嵌入、序列窗口与四项生化特征 | 为比较候选提供互补信息 |
| 实验标注有限，且部分证据只有相对排序 | 冻结预训练编码器，联合学习分类标签与位点排序 | 用小规模任务参数适配现有证据 |
| 需要理解推荐原因并交给后续计算人员 | 输出过滤原因、结构检查、PDB、评分报告及能垒请求 | 让推荐结果能够被复核、复算和继续评估 |
| 后续实验与模拟结果需要持续积累 | 区分湿实验、模拟、审核状态和数据划分 | 为有来源记录的模型迭代提供接口 |

## 核心特色

### 1. 从序列到候选结构与排序的完整设计流程

输入一条标准氨基酸序列后，系统自动枚举新增一个 Cys 的所有位置，对重复序列去重，再调用本地 ESMFold 为每条不同候选预测配套结构。随后提取特征、完成评分与结构检查，输出至多三条不同的推荐序列。

这是单残基插入：长度 N 的输入对应长度 N+1 的候选。系统保留新增 Cys 的位置和原序列映射，确保每个分数对应具体候选、具体位点和对应结构。也支持研究者自行提供候选 PDB，接入已有结构准备流程。

### 2. 序列、结构与生化信息联合建模

模型融合三类信息：CA-only ProteinMPNN 提取的 **128 维结构嵌入**、覆盖目标位点及其前后各 5 个残基的 **231 维序列窗口编码**，以及二硫键状态、溶剂可及表面积（SASA）、前位 Pro 标记和邻近原子数 **4 项生化特征**，形成 363 维位点表示。

结构表征描述空间邻域，序列窗口提供局部残基上下文，生化特征补充暴露程度与邻域环境。项目在预训练结构编码能力之上实现任务适配、候选筛选与结果管理，便于围绕具体位点开展分析。

### 3. 以 364 个可训练参数适配小样本任务

项目冻结 ProteinMPNN 编码器，仅训练一个线性分类头：363 个输入权重加 1 个偏置，共 **364 个可训练参数**。这种方案控制了任务训练的参数规模，也便于复现、特征消融和版本对照。

当前模型同时接收两种监督：已测位点的分类结局，以及 SST 中 Cys3 优于 Cys14 的相对排序。排序监督不要求把较低优先级位点标成阴性；100 个构象共同分摊一条实验观察的权重，保留“多构象、单条实验依据”的数据语义。轻量化指任务分类头，自动结构预测仍需要 ESMFold 的计算资源。

### 4. 推荐结果附带可核查的结构证据

模型评分与显式过滤共同参与筛选。系统记录二硫键与前位 Pro 状态；插入设计还检查结构序列一致性、重原子完整性、关键键长及严重原子重叠，并保存几何证据与过滤原因。

设计流程优先推荐通过过滤的候选。数量不足时，可从具有有效分数的未通过过滤项中补足，同时保留 `eligible=false`、`selection_status=fallback` 及原因提示。缺失结构、序列不匹配或无有效分数的候选不会被赋予虚构分数。评审者可以同时查看推荐清单与全部候选记录，了解候选进入或未进入 Top 3 的依据。

### 5. 网页演示与科研复现使用同一套能力

网页提供序列输入、任务进度、断线重连、Top 3 展示、结构下载、报告与外部评估包下载；命令行覆盖训练、评分、候选设计和反馈处理；Notebook 则串联数据核对、训练、分组诊断及模型结果对照。

ESMFold 在本地完成结构预测，序列不上传到结构预测服务。首次使用需下载权重，缓存完整后支持离线加载。网页同一服务一次处理一个任务，运行日志与结果文件落盘保存，便于演示后继续检查。

### 6. 将候选筛选与后续实验、计算衔接

候选交付包含序列、插入位置、模型分数、配套 PDB、结构来源及能垒计算请求。外部人员可据此准备 GFN2-xTB 或 DFT 计算，再通过反馈校验与独立能垒回归接口接回结果。

湿实验分类反馈与模拟能垒反馈分别管理。只有经人工审核、划入训练集的湿实验反馈可以用于分类模型迭代；模拟数据进入独立接口，电子能垒与活化自由能分别建模。目前已实现交接和训练接口，尚无真实能垒反馈或已验证的能垒模型。

### 7. 将数据语义、版本与复现信息纳入交付

项目使用独立监督掩码区分已测位点与未标注位点，避免把未知结局当作阴性。训练与反馈流程检查重复记录、重复序列位点和留出母蛋白冲突，并提供按来源组/母蛋白组织的诊断与特征消融。

运行输出保留输入与模型哈希、模型版本、筛选策略和完整位点评分；插入设计还记录结构来源与相关代码哈希。正式权重、配置、训练快照和历史对照共同构成可追溯的交付材料。

## 成果与验证

### 已完成的序列设计演示：KTTKS

项目保存了一次真实本地 ESMFold 运行：以 `KTTKS` 为输入，完成全部 **6 个插入位置**的结构预测、特征提取、评分和外部评估包导出。三条推荐序列均通过当次采用的过滤规则。

| 排名 | 推荐序列 | 新增 Cys 位置（1-based） | 模型分数 |
| --- | --- | --- | --- |
| 1 | `KTTCKS` | 4 | 0.9625 |
| 2 | `KTTKCS` | 5 | 0.9455 |
| 3 | `KCTTKS` | 2 | 0.9421 |

可直接查看 [候选清单](results/kttks_esmfold_20260925/results/top3.csv)、[结构报告](results/kttks_esmfold_20260925/results/report.md)和[运行元数据](results/kttks_esmfold_20260925/results/run.json)。这些材料验证了工程流程能够产出可交接的候选；该案例尚未执行 GFN2-xTB 或湿实验验证，分数不能解释为环化效率。

### 当前正式模型：分类与排序联合训练

训练数据包含 **19 条分类记录、22 个监督位点（8 正、14 负）**，以及一条由提供方给出的 SST 排序观察。模型固定训练 50 轮，未根据排序结果追加训练或选择权重。

| 指标 | 结果 | 解释范围 |
| --- | --- | --- |
| 分类训练回代 | 22/22 正确 | 对已有分类监督的拟合结果 |
| SST 训练构象中 Cys3 > Cys14 | 99/100 | 对已纳入训练的排序观察的拟合结果 |
| SST 按轨迹平均后的排序符合数 | 10/10 | 同一案例的构象汇总，非 10 条独立实验 |
| SST 两位点平均分 | Cys3：0.983623；Cys14：0.769741 | 当前训练构象上的优先级比较 |
| 留出整个 SST 来源组后的排序符合数 | 26/100 构象 | 移除相关监督后的分组诊断 |

100 个构象对应同一条排序观察，不能视为 100 个独立实验。当前结果说明模型能够拟合已有分类与排序信息；尚不足以估计新游离短肽的总体准确率或插入设计成功率。详见 [当前训练报告](logs/retrain_sst_20260924/TRAINING_REPORT.md)和 [Model Card](models/MODEL_CARD.md)。

### 工程验证与历史对照

- **自动化验证：** [2026-09-25 自动流程记录](tmp/docs/AUTOMATIC_DESIGN.md)记载 40 项测试通过，并完成短肽 `AG` 的真实结构预测、评分与 ZIP 导出；其几何异常按规则显示警告。该记录反映当时的验证状态。
- **训练与重载复核：** 提供从头执行的 [完整流程 Notebook](notebooks/01_reproducible_workflow.ipynb)，用于训练、分组诊断及正式权重与重训权重的结果对照。[早期验证记录](logs/verification_20260924/verification.json)保存了旧版模型的 28 项测试、CPU 重载与 Notebook 执行结果。
- **保留未改善排序的探索结果：** 旧版模型在还原封端 SST 的 100 个 MD 生产快照中仅 7 个给出 Cys3 更高；固定 λ=0.5 的 PROPKA 修正为 6/100，未改善排序。这些记录作为方法迭代依据保留，不能与当前模型的训练回代混为独立验证。

已知适用边界包括：训练规模小、缺少独立阴性测试、历史结构特征计算材料不完整；自动插入设计目前每个候选使用一个未松弛预测构象，尚未加入 MD 多构象聚合或反应建模；端基修饰未进入序列编码，前位 Pro 过滤仍是待核验的保守假设。后续验证应冻结模型，用新的来源组、预先确定的候选与对照开展统一计算和实验。

## 技术方案

### 流程架构

```text
已有 Cys 位点评估                         新多肽的单 Cys 插入设计
序列 + Cα 坐标 + 生化特征                 输入序列 → 枚举 N+1 个插入操作
             │                                       │
             │                           序列去重 → 本地 ESMFold / 配套 PDB
             │                                       │
             │                               结构校验与生化特征提取
             └───────────────────┬───────────────────┘
                                 │
            ProteinMPNN 结构嵌入 128 维（冻结编码器）
              + 局部序列窗口 231 维 + 生化特征 4 维
                                 │
                        363 维融合 → 线性头 → Sigmoid
                                 │
                   按任务执行阈值、过滤、排序与去重
                                 │
               Top 3 + 全量评分 + 结构证据 + 运行元数据
                                 │
               外部实验 / 能垒计算 → 审核归档 → 独立迭代接口
```

### 模型训练与监督方式

结构编码器采用 CA-only ProteinMPNN。11 残基序列窗口使用 21 类编码，边界以 `X` 填充；四项生化特征按固定尺度归一化。标签、样本来源与实验结局不作为预测特征。

当前训练的联合损失为：

```text
L = [Σ BCE（22 个分类位点）
     + mean（100 个构象上的 softplus(-(z_Cys3 - z_Cys14))）] / 23
```

其中 `z` 为分类头输出的 logit。每个分类位点权重为 1，整条 SST 排序观察的权重为 1；100 个构象只作为该观察的输入增广。未标注 Cys 不参与分类损失，Cys14 不因排序较低而被改标为阴性。

| 配置项 | 默认设置 |
| --- | --- |
| 正式模型版本 | `sst-ranking-20260924-v1` |
| 正式权重 | `models/site_predictor.pt` |
| 结构编码器 | 冻结的 CA-only ProteinMPNN |
| 分类头 | 363 维输入，364 个可训练参数 |
| 优化器 | AdamW |
| 学习率 / 权重衰减 | 0.01 / 0.01 |
| 训练轮数 / 随机种子 | 50 / 20260916 |
| 分类判定阈值 | 0.5 |
| 排序监督配置 | `data/ranking/sst_ordering.json` |

配置见 [configs/train.json](configs/train.json)。最终模型采用固定训练轮次，不依据留出诊断调整阈值或选择权重。移除配置中的 `ranking` 键可复现纯分类训练，需另用新的输出目录。

### 两类任务的筛选规则

| 规则 | 已有位点评分：`predict.py` / `screen.py` | 插入设计：`design.py` / 网页 |
| --- | --- | --- |
| 评分对象 | `feature_mask=1` 的已有 Cys | 每个插入操作中的新增 Cys |
| 0.5 阈值 | 通过过滤且分数 ≥ 0.5 才入选 | 单独记录分类判定，不作为设计入选硬门槛 |
| 主要过滤 | 二硫键状态、前位 Pro | 结合既有规则、N 端限制与结构几何检查 |
| 排序 | 分数降序；同分按记录编号和位点位置 | 通过过滤者优先，组内分数降序；同分按候选编号 |
| Top 3 不足 | 输出实际合格数量 | 可用有效评分的过滤项补足，明确标记原因 |
| 去重 | 按位点输出 | 按完整候选序列去重，优先保留通过过滤的位点解释 |

SASA 与邻近原子数参与模型输入，不设置额外硬阈值。插入设计即使启用补足，也可能因不同序列或有效分数不足而少于三条。

### 插入候选的结构测量

- **SASA：** 使用 Biopython Shrake–Rupley，探针半径 1.4 Å、960 个采样点，取目标 Cys 的 CB 与 SG 重原子面积之和。
- **二硫键几何标记：** 目标 SG 与其他 Cys SG 距离小于 2.5 Å；距离标记本身不是化学成键状态的独立证明。
- **邻近原子数：** 统计目标 SG 周围 5 Å 内的非氢原子，仅排除该 SG 本身。
- **结构质量检查：** 检查主链与目标 Cys 的关键键长，以及非相邻残基间小于 1.2 Å 的严重原子重叠。

手工输入要求候选序列与 PDB 完全对应，使用分离单链、标准氨基酸、无配体或水分子且重原子完整的结构。原始序列的 PDB 不能替代插入后结构。详见 [插入设计说明](tmp/docs/DESIGN_WORKFLOW.md)。

## 快速开始与评审演示

建议先查看上述 KTTKS 候选和结构报告，再运行已有位点评分；安装结构预测依赖后，可进一步体验网页设计流程。

### 1. 配置环境

仓库记录的已验证环境如下，完整依赖见 [requirements.txt](requirements.txt)。

| 项目 | 配置 |
| --- | --- |
| 操作系统 | Windows 11 x64，build 26200 |
| Python | 3.13.5 |
| PyTorch / NumPy | 2.11.0+cu128 / 2.5.2 |
| GPU | NVIDIA GeForce RTX 4050 Laptop GPU，6 GB |
| NVIDIA 驱动 / PyTorch CUDA 运行时 | 592.82 / 12.8 |

在仓库根目录创建独立环境：

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements.txt
```

下文 `python` 指已安装项目依赖的解释器，例如 Windows 下的 `.venv/Scripts/python.exe`。位点评分使用 CPU；训练的结构编码阶段自动选择可用 CUDA 或 CPU，分类头在 CPU 上优化。其他平台需使用对应的 PyTorch 分发包并重新核验，Linux 尚未实机验证。

### 2. 运行已有位点评分

```bash
python predict.py --input data/examples/input.jsonl --output results/review.csv --track "实际参赛赛道"
```

| 生成文件 | 内容 |
| --- | --- |
| `results/review.csv` | Top 3 位点、分数、结构特征及版本信息 |
| `results/review.sites.jsonl` | 全部已评分位点、分类判定与过滤原因 |
| `results/review.run.json` | 输入/模型哈希、软件版本、阈值及候选数量 |

可先查看随附的 [示例候选清单](results/results.csv)。示例输入来自训练数据，用于运行演示，不是独立测试集。

### 3. 通过网页设计新序列

**当前目录兼容性提示：** 专题文档已迁移至 `tmp/docs/`，但 `src/folding.py` 仍从 `docs/BARRIER_INTERFACE.md` 复制交接说明，自动流程会在评估包导出阶段失败。运行完整自动设计前需同步该路径；下述操作说明及已保存演示反映已实现的功能。

```powershell
python -m pip install -r requirements-folding.txt
python app.py
```

保持服务终端运行，在浏览器打开 [本地操作界面](http://127.0.0.1:8765)。

1. 在“目标多肽序列”中输入如 `AGKSTV` 的序列。仅接受 20 种标准氨基酸单字母代码，不含 FASTA 标题、编号或端基修饰；自动忽略空白并转为大写，默认长度为 2–100 个残基。
2. 点击“开始预测”，通过“查看运行详情”观察进度。同一服务一次处理一个任务；连接中断后可点击“重新连接任务”，服务仍在运行时也可刷新恢复。
3. 在“推荐序列 / TOP 3”中查看高亮的新增 `C`、候选分数和提示。新增位置按候选序列从 1 开始编号，候选长度比输入多 1。
4. 下载候选 PDB、Markdown 报告或外部评估 ZIP，供后续分析与计算交接。

首次使用需下载 ESMFold 权重，缓存完整后可执行 `python app.py --local-files-only`。默认 `--device auto` 仅在可用 GPU 显存至少 12 GiB 时选择 CUDA，否则使用 CPU；CPU 仍需要足够内存，自动结构预测可能耗时较长。此资源要求与轻量位点评分不同。

结果默认保存在 `results/web_jobs/<任务编号>/design/`，日志位于同一任务目录的 `progress.log`。服务重启后文件保留，但网页不能恢复旧任务；可直接读取磁盘结果。若预测未完成，应检查依赖、权重下载或内存相关日志后重新提交。完整参数见 [自动结构预测与多肽输入窗口](tmp/docs/AUTOMATIC_DESIGN.md)。

### 4. 命令行设计与手工结构接入

自动结构预测与筛选使用同一入口：

```powershell
python design.py run --sequence AGKSTV --record-id demo_only --output results/auto_design_demo
```

也可分步枚举候选，再提供每条插入后序列的真实配套 PDB：

```powershell
python design.py prepare --sequence AGKSTV --record-id demo_only --output results/design_demo_requests
# 将生成的结构清单模板填写并另存为 structures.json
python design.py score --candidates results/design_demo_requests/candidates.jsonl --structures results/design_demo_requests/structures.json --output results/design_demo_scored
```

`AGKSTV` 仅为流程示例。默认要求所有候选都有结构；仅评估部分候选时需显式使用 `--allow-partial`，输出会记录覆盖范围。

### 5. 复现训练与完整流程

```bash
python run.py --output logs/demo --track "实际参赛赛道"
```

该命令完成模型训练与参数优化、位点评分和候选筛选，保存权重、损失日志、分组评估、候选清单及运行耗时。默认训练 50 轮，可加 `--epochs 1` 检查环境；此入口不执行 ESMFold 自动插入设计。

**路径约定：** 配置及命令行中的运行相对路径以仓库根目录为基准；手工结构清单中的 PDB 相对路径以清单所在目录为基准。输出文件或训练目录须使用新路径，避免覆盖已有结果。默认预测文件 `results/results.csv` 已存在，复跑时应另设输出；`--track` 填写实际参赛赛道，未填写则标记为演示。

## 训练复现与反馈扩展

### Notebook 与独立入口

[完整流程 Notebook](notebooks/01_reproducible_workflow.ipynb)包含中文步骤说明、数据核对、50 轮训练、分组诊断、正式模型与重训模型的 Top 3 对照及原始评分一致性检查，并保存了执行输出。

```bash
python -m pip install -r requirements.txt
python -m jupyterlab
```

选择相同项目环境的 Python 内核，执行 **Restart Kernel and Run All Cells**。默认 50 轮，参数单元格可改为 1 轮冒烟检查。可从仓库根目录或 `notebooks/` 启动；内核工作目录位于仓库外时，先设置 `CYSRANK_ROOT`。每次执行写入新的 `logs/notebook/<运行编号>/`，不覆盖正式权重。Notebook 的候选同样属于训练回代演示。

独立训练与筛选命令：

```bash
python train.py --config configs/train.json --output logs/reproduction
python screen.py --input data/examples/input.jsonl --output results/screen.csv --track "实际参赛赛道" --top-k 3
```

`screen.py` 与 `predict.py` 使用相同规则，候选不足时输出实际合格数量。工程测试入口：

```bash
python -m unittest discover -s tests -v
```

### 湿实验反馈

```bash
python feedback.py --input data/feedback/incoming.jsonl --output data/feedback/reviewed.jsonl
python train.py --feedback data/feedback/reviewed.jsonl --output logs/feedback_iteration
```

反馈文件由使用者提供。校验与归档保留观察来源、审核状态和划分信息；只有已审核且属于训练集的湿实验记录进入分类训练，模拟、留出与未审核记录独立保留。字段要求见 [反馈接口说明](tmp/docs/FEEDBACK.md)。

### 外部能垒反馈

`barrier.py requests` 导出请求，外部计算人员补充结果后使用 `feedback.py` 校验，再由 `barrier.py train` 训练独立回归器、`barrier.py predict` 预测同一计算协议下的新位点。支持 GFN2-xTB 或 DFT 反馈，计算失败和未审核结果不参与训练。

网页和设计入口不执行量化能垒计算。交接 PDB 还需由计算方定义试剂、反应步骤、端基/质子化、电荷、溶剂、参考态与过渡态核验；未计算能垒保持为空，不以结构优化能或反应能代替。详见 [能垒接口说明](tmp/docs/BARRIER_INTERFACE.md)。

## 数据与结果格式

### 已有位点评分输入

输入为 UTF-8 JSONL，每行一条序列。设长度为 N，评分输入至少含 3 个残基。

| 字段 | 形式 | 含义 |
| --- | --- | --- |
| `record_id` | 唯一非空字符串 | 输入记录编号 |
| `seq` | 长度 N 的字符串 | 氨基酸序列 |
| `xyz` | N × 3 数组 | Cα 坐标，单位 Å |
| `is_disulfide` | N 个 0/1 | 二硫键状态 |
| `sasa` | N 个非负数值 | 溶剂可及表面积，单位 Å² |
| `prev_is_pro` | N 个 0/1 | 前一残基是否为 Pro |
| `nearby_atoms` | N 个非负数值 | 5 Å 范围内非氢原子数 |
| `feature_mask` | N 个 0/1 | 结构特征已准备、需要评分的位点 |
| `structure_file` | 可选相对路径 | 仓库内对应结构文件 |

模型只评分 `feature_mask=1` 的 Cys，推理无需标签。SASA 为 0 可以是有效测量，缺失特征需补齐，不能以 0 代替。结构文件应与记录真实对应。数据来源与划分见 [数据说明](data/README.md)。

### 标准输出与交付包

已有位点的候选 CSV 包含候选编号、赛道、序列、Cys 位置、分数、排名、模型及运行版本、权重哈希、结构引用、四项生化特征和备注。位点采用 1-based 编号，标准评分 CSV 的分数保留 8 位小数。

插入设计进一步输出 `top3.csv`、`all_candidates.jsonl`、`report.md`、标准化结构副本、实际评分的模型输入、能垒请求及 `run.json`。自动流程的外部评估 ZIP 汇集候选、结构、来源与交接文件，便于交给后续计算人员。

## 项目结构

```text
CysRank-GCLSC/
├── README.md                 # 项目介绍、评审演示与技术说明
├── requirements.txt          # 基础运行与 Notebook 依赖
├── requirements-folding.txt  # 本地 ESMFold 依赖
├── requirements-md.txt       # 历史 MD 案例依赖
├── requirements-propka.txt   # PROPKA 对照依赖
├── configs/                  # 分类与排序联合训练配置
├── data/                     # 数据、监督掩码、排序观察与来源记录
├── src/                      # 特征、训练、推理、结构与设计逻辑
├── models/                   # 正式权重、预训练权重及 Model Card
├── app.py                    # 本地网页服务
├── web/                      # 序列输入与结果展示界面
├── train.py                  # 训练入口
├── predict.py                # 标准化位点评分入口
├── screen.py                 # 候选筛选入口
├── run.py                    # 训练到筛选的完整流程
├── design.py                 # 单 Cys 插入、结构预测与评分
├── barrier.py                # 能垒请求、独立训练与预测
├── feedback.py               # 反馈校验与归档
├── results/                  # 候选结果、结构报告与历史对照
├── logs/                     # 训练、分组评估与验证记录
├── notebooks/                # 可执行演示与案例分析
├── tests/                    # 工程测试
├── tools/                    # 案例复算与本地打包工具
├── tmp/docs/                 # 本地专题资料，当前不纳入提交包
└── licenses/                 # 第三方许可证
```

## 专题资料与交付说明

### 案例复算与历史记录

还原态 SST 案例按提供方确认的 Ac/NH2 端基、还原巯基、pH 7.5 和 300 K 完成 10 起点短时显式水模拟。旧版模型平均分为 Cys3=0.8918、Cys14=0.9853，按起点平均均为 Cys14 更高，未稳定复现实验排序。相关脚本绑定旧版权重，以保留版本对照。

可查看 [案例说明](tmp/docs/REDUCED_PEPTIDE_CASE.md)、[案例 Notebook](notebooks/02_reduced_peptide_case.ipynb)、[逐位点评分](results/2mi1_reduced_capped_ph75_20260924/analysis/results.csv)和 [PROPKA 对照报告](results/2mi1_propka_lambda05_20260924/REPORT.md)。仅复算随附结构可执行：

```bash
python tools/analyze_2mi1_reduced_md.py --output results/2mi1_reanalysis
```

重做模拟需另装 `requirements-md.txt`，OpenCL 设置见案例说明。该流程不插入残基、不计算反应能垒，也不改变正式模型。[KTTKS 文献案例](tmp/docs/KTTKS_CASE_STUDY.md)和 `results/kttks_design_requests/` 保存前期候选准备材料；已完成的自动预测结果位于 `results/kttks_esmfold_20260925/`。

### 本地打包与目录兼容性

专题文档现位于 `tmp/docs/`，历史版本、重复运行产物与过程文件归档至 `tmp/`。[历史整理说明](tmp/docs/PROJECT_LAYOUT.md)中的根目录 `docs/` 描述反映迁移前布局。

打包工具按明确的文件与目录范围生成 ZIP 和逐文件 SHA-256 清单，排除 `tmp/`、`.git/`、本机环境、缓存及重复 Notebook 运行目录，只生成本地包、不上传。**当前 `tools/package_submission.py` 仍要求根目录存在 `docs/`，需同步目录配置后再运行；同时应确认评审所需专题文档已纳入交付范围。**

```bash
python tools/package_submission.py --output tmp/agents/submission/CysRank-GCLSC.zip
```

输出须为新路径。历史训练到筛选流程耗时约 12–37 秒，不包含 ESMFold 结构预测；实际耗时随硬件与加载状态变化。

### 数据来源与第三方声明

CysRank 基于 ProteinMPNN 的预训练结构编码能力，实现 Cys 位点任务适配、生化特征融合、候选设计筛选与反馈管理；自动结构预测使用 ESMFold。ProteinMPNN 的 MIT 许可保留于 [licenses/ProteinMPNN-MIT.txt](licenses/ProteinMPNN-MIT.txt)，本地版本与哈希见 [第三方说明](tmp/docs/THIRD_PARTY.md)，结构预测相关依赖说明见 [自动流程文档](tmp/docs/AUTOMATIC_DESIGN.md)。

部分原始实验记录、数据授权、获取时间及结构特征复算材料尚需提供方补充。正式提交范围与候选字段应核对实际赛道要求，详见 [数据说明](data/README.md)和[提交要求核对表](tmp/docs/SUBMISSION_CHECKLIST.md)。
