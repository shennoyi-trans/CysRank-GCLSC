# 结果说明

results.csv 为随附模型对 data/examples/input.jsonl 的训练回代演示，最多保留通过 0.5 阈值和硬过滤的全局 Top 3。赛道“未指定（演示）”需正式填写，不得冒充独立验证结果。

candidate_id、track、sequence、prediction_score、model_version、run_version、remarks 对应候选编号、赛道、序列、关键指标、模型/运行版本和备注。另提供位点、结构引用、结构特征和模型哈希。分数保留 8 位小数，缺失全原子结构时 structure_file 为空并在备注说明；SASA 单位 Å²，位点为 1-based。

.sites.jsonl 保留全部打分、硬过滤和拒绝原因；.run.json 记录输入及模型哈希、环境版本、阈值、排序和候选数量。不提供未经验证的置信区间。

kttks_design_requests/ 保留 KTTKS 的六条单 Cys 插入候选和结构交接模板，目前无结构、模型评分或 Top3。案例依据见 ../docs/KTTKS_CASE_STUDY.md。


## 还原态测试肽案例

2mi1_reduced_capped_ph75_20260924/ 保存本轮成功的结构生成与评分结果。analysis/results.csv 是 100 个生产快照、每帧两个已有位点的 200 行清单，包含实际结构文件路径；analysis/comparison.csv 还包含最小化和平衡结构。它不是 200 个新设计候选，也不是独立测试。

实验期望 Cys3 > Cys14，当前仅 7/100 帧符合，10 个起点的平均排序均相反。全部快照保留，不挑选有利构象。参数与适用范围见 [正式案例说明](../docs/REDUCED_PEPTIDE_CASE.md)，系统与能量等原始模拟日志在 logs/reduced_md_20260924。失败、中断及仅改标记的早期尝试已归档至 tmp/archive/reduced_case_process_20260924。


## 当前模型：SST 排序监督版

正式权重已更新为 sst-ranking-20260924-v1。原 22 个分类标签不变，新增一条 Cys3 > Cys14 实验排序监督；100 个构象共同占一条观察的权重。SST 训练回代 99/100 帧排序符合，不能作为独立验证。此前还原态/PROPKA 报告保留为旧模型对照。见 [重训说明](../docs/SST_TRAINING.md)。
