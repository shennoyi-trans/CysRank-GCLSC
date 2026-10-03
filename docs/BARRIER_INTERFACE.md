# 外部能垒计算与训练接口

计算人员负责量化计算；本项目负责导出计算请求、校验结果、训练与预测。接口版本为 `cysrank-barrier-v1`，统一交换 UTF-8 JSONL，每行一个序列位点的已汇总模拟结果，不解析特定版本的 xTB/ORCA 日志。

## 交接流程

```powershell
# 使用已有结构特征输入导出待计算请求，文件不能已存在
python barrier.py requests --input data/examples/input.jsonl --output tmp/agents/barrier_requests.jsonl

# 计算方填写请求中的空字段，完成审核后交回 external_results.jsonl
python feedback.py --input external_results.jsonl --output data/feedback/barriers_round1.jsonl

# 单独训练模拟能垒回归器，不改变已有实验分类头
python barrier.py train --input data/feedback/barriers_round1.jsonl --output logs/barrier_round1

# 使用本轮生成的计算协议预测新候选；所有路径相对仓库根目录
python barrier.py predict --input candidates.jsonl --model logs/barrier_round1/barrier_model.json --protocol logs/barrier_round1/protocol.json --output results/barrier_predictions.jsonl
```

示例输入来自训练数据，仅展示接口操作。请求不含虚构能垒；`pending`、空元数据请求不能进入反馈归档或训练。正式使用时应换成实际候选及其真实结构特征。当前接口不执行 Cys 插入、结构生成或量化计算。

## 计算方填写什么

请求自带 `sample`（seq、Cα 坐标、四项生化特征）、1-based Cys 位置和 `sample_sha256`。这些描述模型输入的反应前结构，不应替换成含产物或过渡态信息的特征。计算方若改变模型输入结构，需重新导出请求；实际全原子计算结构另外通过 `structure_reference` 关联。

公共字段沿用 `feedback.py` 接口：

| 字段 | 约定 |
|---|---|
| observation_id | 一次交付的唯一 ID，重复观察需人工合并 |
| evidence_type / outcome | 必须为 `simulation` / `null`；不把低能垒映射成实验阳性 |
| parent_group | 同一母蛋白/同源家族使用同一组标识，插入变体也应归在原组 |
| split | `train` 或 `holdout`，必须在训练前确定 |
| assay_id / conditions | 计算批次及非空实验/模拟条件对象 |
| provenance | 原始日志或可追溯报告路径/URI |
| reviewed | 人工审核后的布尔值 |

新增的 `barrier` 对象：

| 字段 | 约定 |
|---|---|
| schema_version | `cysrank-barrier-v1` |
| status | `completed` 或 `failed`；请求初始为 `pending` |
| value / unit | 已计算的能垒差值；支持 `kcal/mol`、`kJ/mol`、`eV`、`Hartree`，内部统一为 kcal/mol |
| quantity | `delta_e_dagger` 电子能垒，或 `delta_g_dagger` 活化自由能；不是反应能或总能量 |
| quality_passed / quality_notes | 是否通过计算方质量检查及依据；布尔值必须真实填写 |
| structure_reference | 实际全原子计算结构、原子映射或结构清单路径/URI |
| sample_sha256 | 导出请求时自动生成，不需计算人员手写 |
| charge / multiplicity | 计算体系总电荷（整数）及自旋多重度（正整数） |
| protocol | 以下统一计算方案 |
| failure_reason | 失败时必填；失败必须 `value=null, quality_passed=false` |

`protocol` 的必填项：

- `software`、`software_version`、`method`：例如 xTB、实际版本、GFN2-xTB；也支持 DFT 方法名称。
- `reaction_id`：确定 CA 衍生物和反应定义的稳定标识。
- `reaction_step`：例如指定的加成步或成环步；多步不能任意取最大值后称为总体能垒。
- `reference_state`、`standard_state`：例如相对哪个反应物/复合物，以及标准态处理约定。
- `solvent`、`temperature_kelvin`：包含溶剂模型及溶剂组成说明，温度为正数 Kelvin。
- `protonation_policy`、`conformer_policy`：质子化与构象采样/汇总规则。
- `settings_reference`：完整计算方案的稳定版本或哈希；其中应说明反应路径、过渡态与频率/IRC 检查、收敛设置、热修正等。不能填每个样本都不同的日志路径；样本日志填 provenance。

计算人员负责从软件输出中得到约定差值，并确认搜索到了对应反应路径。质量检查失败、未收敛、只得到未经核验的路径最高点时应设 `quality_passed=false`。接口只能检查数据格式和标记，不能替代化学审核、验证外部日志真实性或自动判断过渡态。

单位转换不会混淆能垒定义。有限的负值允许保留以供审核，不自动截断为零；异常负值应由计算方检查参考态和修正方案。

## 训练与模型边界

- 仅 `reviewed=true`、`completed`、`quality_passed=true` 且 `split=train` 的数据训练。
- 留出记录只用于 MAE/RMSE；未提供留出集时报告明确标明，不声称泛化有效。
- 一个模型只能学习一种 quantity 和完全一致的 protocol。方法、反应步骤、溶剂或参考态不同必须分开训练。协议 JSON 中额外字段也参与一致性检查。
- 同一序列/位点/协议只允许一条通过审核的汇总结果，不把重复计算或多个构象当作多个独立样本。如何汇总由计算方预先在 conformer_policy 中定义。
- 拒绝训练/留出母蛋白组交叉和序列包含关系交叉；完整同源性和变体归组仍需人工审核。
- 回归器复用冻结 ProteinMPNN 的 363 维特征，采用带 L2 正则的线性回归（ridge）。特征缩放只从训练集拟合，默认 alpha=1，不根据留出集自动调参。
- 第一版至少要求两个不同训练位点以使接口可运行，这不是足以证明预测质量的样本量标准。输入暂不显式编码 CA、pH 或总电荷；只能在固定反应协议和一致结构准备政策下使用，不支持任意条件迁移。
- 输出是 `predicted_barrier_kcal_mol`，按数值升序排列，不是概率、产率或实验标签；不自动融合到旧分类分数，也不自动更改硬过滤规则。
- 每轮使用同一协议的累计、去重反馈重新拟合；不同协议分目录管理。保存模型、协议、输入哈希、接受的记录、预测及误差报告。已有模型和输出不覆盖。

## 交付给计算人员的简要约定

我们提供序列、目标 Cys、模型输入结构特征及其哈希。计算方准备实际量化体系并交回：能垒差值、单位、能垒种类、完整协议、结构/日志来源、成败状态和质量检查结论。计算失败返回失败记录而不是零值或阴性标签。统一协议尚未确定时可以先填写请求模板，但不能用占位符训练正式模型。
