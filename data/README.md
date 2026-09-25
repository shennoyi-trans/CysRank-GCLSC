# 数据组织与版本管理

2026-09-24 导入来源方“第二部分任务 (2)”；原始来稿保留在 tmp 中。正式 raw/ 含 8 条正例与 11 条负例，22 个监督位点（8 正、14 负）。沿用既定 Pro 冲突隔离，排除 VNRRPCFSALT，不改标阴性；provenance/positive_before_exclusion.jsonl 保存本批次排除前 9 条正例。

原 ca_label 按已测试位点掩码解释，正负结局由文件来源决定；转换后保存 tested_mask、label_mask、confirmed_negative_mask。3 个未标注 Cys 不参加损失。实验来源沿用用户说明，原始实验记录未独立审计。

| 路径 | 用途 |
| --- | --- |
| raw/ | 本轮训练源文件 |
| processed/training.jsonl | 本轮训练实际保存的监督快照，不重复训练 |
| examples/input.jsonl | 相同记录移除标签后的训练回代示例 |
| external/validation.jsonl | 248 条全零标签的未验证候选，不视为确认阴性，不训练或调参 |
| structures/test.pdb | 沿用 2MI1；还原结构已用于 SST 排序训练，不是独立测试 |
| provenance/parent_sequences.jsonl | 历史完整序列，仅用于精确包含关系分组 |
| provenance/revision_20260924.json | 来源哈希、排除原因和变更范围 |

## 特征更新与评估划分

本轮整合正例 SASA 更新，并保留负例数据以延续既有实验监督。更新范围为片段结构特征，未进行游离短肽构象重建；特征计算的完整复算材料仍待补齐。

按历史母序列精确包含关系组织 5 个分组，配套分组诊断和特征消融。含全部阴性的组留出后训练单类，该折跳过；其余可评估记录用于正例诊断。训练包含测试肽 AGCKNFFW、KTFTSC 片段，2MI1 因此不列为独立测试。版本间样本组成不同，不将指标差异归因于单项特征更新。

原始数据、训练快照与源文件哈希可供追溯；实验来源依据提供方说明。来源许可、完整结构计算记录和同源/预训练重叠核验尚需补齐。旧数据与模型保存在 tmp/archive/before_part2_retrain_20260924。


## 当前模型：SST 排序监督版

正式权重已更新为 sst-ranking-20260924-v1。原 22 个分类标签不变，新增一条 Cys3 > Cys14 实验排序监督；100 个构象共同占一条观察的权重。SST 训练回代 99/100 帧排序符合，不能作为独立验证。此前还原态/PROPKA 报告保留为旧模型对照。见 [重训说明](../docs/SST_TRAINING.md)。
