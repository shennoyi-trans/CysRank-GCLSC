# 结果说明

results.csv 为随附模型对 data/examples/input.jsonl 的训练回代演示，最多保留通过 0.5 阈值和硬过滤的全局 Top 3。赛道“未指定（演示）”需正式填写，不得冒充独立验证结果。

candidate_id、track、sequence、prediction_score、model_version、run_version、remarks 对应候选编号、赛道、序列、关键指标、模型/运行版本和备注。另提供位点、结构引用、结构特征和模型哈希。分数保留 8 位小数，缺失全原子结构时 structure_file 为空并在备注说明；SASA 单位 Å²，位点为 1-based。

.sites.jsonl 保留全部打分、硬过滤和拒绝原因；.run.json 记录输入及模型哈希、环境版本、阈值、排序和候选数量。不提供未经验证的置信区间。
