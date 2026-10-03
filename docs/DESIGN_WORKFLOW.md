# 陌生多肽的 Cys 插入设计

本流程实现真正的单残基插入：输入长度 N，候选长度 N+1，删除新增的 C 后必须还原原序列。纯序列输入生成待评估候选；提供每条插入后候选的配套全原子结构后，执行评分、结构检查、Top3 导出及能垒计算请求交接。

原始多肽的 PDB 不等于插入后候选的 PDB。现可使用 `python app.py` 的序列输入窗口或 `design.py run` 自动调用本地 ESMFold，为每条不同候选预测配套结构；见 [自动流程](AUTOMATIC_DESIGN.md)。以下 prepare/score 接口仍支持手工提供结构。结构必须与候选序列完全一致。

## 第一步 输入序列

安装 `requirements.txt` 中的依赖。新增加的 Biopython 1.88 用于 PDB 读取和 SASA 计算。

```powershell
python design.py prepare --sequence AGKSTV --record-id demo_only --output results/design_demo_requests
```

`AGKSTV` 是流程说明用序列，不是已验证的设计案例。正式使用时替换为真实陌生多肽。接受大写的 20 种标准氨基酸单字母代码，至少 2 个残基。

上述功能测试的历史结果已归档至 tmp；正式目录保留 results/kttks_design_requests 下的文献案例。再次运行时请选择新的输出目录。

输出：

- `candidates.jsonl`：全部 N+1 个插入操作、新序列、新 Cys 位置和原序列映射。
- `candidate_sequences.fasta`：按序列去重的结构预测交接文件。
- `structure_manifest.template.json`：待填写的候选结构清单。
- `preparation.json`：本次生成范围与待评估状态。

插入位置采用“原序列第几位之后”：0 表示 N 端，N 表示 C 端。新增 C 的候选序列编号为该数值加 1。相邻已有 C 的不同插入操作可能得到相同序列，仍保留独立的位点映射，但只需要一份相同序列的结构。

## 第二步 提供候选结构

将模板另存为 `structures.json`，为每条不同序列填写：

```json
[
  {
    "sequence": "ACGKSTV",
    "structure_id": "demo_only_insert_C_after_1",
    "candidate_ids": ["demo_only_insert_C_after_1"],
    "structure_file": "pdb/demo_only_insert_C_after_1.pdb",
    "chain_id": "A",
    "model_index": 1,
    "structure_method": "填写实际工具、版本及结构准备方法",
    "structure_notes": "填写结构状态、质量信息和来源"
  }
]
```

这只是单条格式示例，未提供虚构 PDB。通常应填写模板中的全部序列。

- PDB 路径相对 `structures.json` 所在目录解析，也接受绝对路径。
- 仅支持分离的单链、标准氨基酸、无配体/水分子的全重原子 PDB；氢原子在特征测量时去除。
- 必须明确链号；多模型 PDB 用 1-based `model_index` 选择构象。不同构象需要分别运行，不会自动挑最有利构象。
- 严格校验序列与每个残基的重原子完整性；不接受未经处理的 alternate locations。
- N 端插入会被生成，但当前反应要求相邻主链酰胺，游离 N 端新增 Cys 不属于通过过滤的候选。合格候选不足时可按分数补入 Top3，并明确提示原因。若实际使用 N 端化学修饰，需要另行扩展化学结构支持。

## 第三步 评分与导出

```powershell
python design.py score --candidates results/design_demo_requests/candidates.jsonl --structures results/design_demo_requests/structures.json --output results/design_demo_scored
```

默认要求所有候选序列都有结构；确实只评估部分结构时可显式加 `--allow-partial`，报告会标明仅覆盖已提供结构，不能宣称得到全部位置中的最优解。每次输出目录须为新目录，避免覆盖历史结果。

输出：

| 文件 | 内容 |
|---|---|
| top3.csv | 至多三条不同序列、插入位置、原始分数和结构文件路径 |
| all_candidates.jsonl | 全候选分数、过滤原因、原子距离和几何证据 |
| report.md | 可阅读的 Top3 结构报告及局限 |
| structures/ | 所选单链构象的标准化重原子结构副本 |
| model_inputs.jsonl | 所有实际已评分位点的模型输入 |
| top3.model_inputs.jsonl | 选中候选的特征输入 |
| barrier_requests.jsonl | 有推荐候选时自动生成，能垒字段为空，供计算人员填写 |
| run.json | 输入/模型/代码哈希、筛选策略、覆盖范围和限制 |

默认 Top3；可用 `--top-k` 调整数量，输出文件名保持 `top3.*` 以保持接口固定。

## 排序与结构证据的含义

新候选复用当前正式位点模型，只对新增 Cys 评分。先执行几何检查和该检查点已有的二硫键/前位 Pro 过滤。通过过滤的候选优先，数量不足时按内部评分补入未通过过滤的候选；各组内按原始分数降序排列，最终按完整序列去重。重复序列优先保留通过过滤的新增位点解释。补入项保留 eligible=false，并在 warning 中提示成功可能性较低及原因。

这里的设计排序不以 0.5 为硬门槛；`classified_positive_at_0_5` 单列保留原模型判定。即使某候选分数低于 0.5，也可能属于当前候选中的相对高分项。过滤候选会用于补足，但不同序列本身不足三条或缺少真实分数时，仅输出实际可用数量。

结构测量口径固定：

- SASA：Biopython Shrake–Rupley，探针半径 1.4 Å，960 采样点，重原子结构上目标 Cys 的 CB+SG 面积和。
- 二硫键几何标记：目标 SG 与其他 Cys SG 距离小于 2.5 Å。距离标记不是对实际化学成键状态的独立证明。
- 邻近原子数：以目标 SG 为中心，5 Å 内非氢原子数，仅排除该 SG 本身。
- 几何检查：所有主链 N–CA、CA–C、C–O、相邻残基 C–N 距离，目标 Cys 的 CA–CB、CB–SG 距离，以及间隔超过一个残基的原子间小于 1.2 Å 的严重重叠。
- 不把 PDB B-factor 自动解释成 pLDDT，也不输出未计算的能量或置信区间。

这些是结构检查和模型优先级，不是实验有效性证明。历史训练 SASA/邻近原子数的原始计算过程不完整，新测量可能与训练分布不同。现有模型没有接受插入设计的独立验证。前位 Pro 硬过滤仍按交付检查点保留，是待核验假设；没有因本次新接口而更改旧模型规则。

## 与能垒闭环衔接

计算人员按 `BARRIER_INTERFACE.md` 完成 `barrier_requests.jsonl` 中的计算设置、实际能垒和质量审核。插入变体自动使用同一 parent_group，防止把同源候选拆入训练与留出；与其他数据合并时仍需人工统一母蛋白/同源家族标识。

当前只导出请求，不计算能垒，也不自动把能垒回归模型融合进位点分类分数。请求所附结构引用由运行机器生成；迁移给其他人时，应随附 `structures/` 并调整引用路径。

## 验证

```powershell
python -m unittest discover -s tests -v
```

测试覆盖真实插入与位置映射、重复序列去重、缺原子与错序列拒绝、二硫键距离、几何异常、完整评分及能垒请求输出。端到端测试结构为临时合成夹具，只用于工程验证，不进入正式训练或作为科学结果保存。

补充：2MI1 的还原态松弛是已有 14 肽的特定案例，不是任意插入候选的自动结构预测入口。见 [案例说明](REDUCED_PEPTIDE_CASE.md)。
