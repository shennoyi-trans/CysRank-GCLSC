# 附件4与反馈意见核对

2026-09-24 编写，2026-10-03 更新依赖与文档位置。依据 tmp/docs/附件4-代码提交要求0721-V2.docx 和反馈意见（3）。正式包提供可运行入口、模型、依赖、数据、日志和案例，不以目录整理替代科学效果验证。

| 要求 | 正式位置与状态 |
| --- | --- |
| 源码、环境、命令 | src、tools；唯一完整 requirements.txt；README 与 docs 使用指南 |
| 训练、种子、最终权重 | train.py、configs、logs/retrain_sst_20260924、models；logs/retrain_20260924 保留旧版对照 |
| 训练/优化/推理示例 | run.py 和 notebooks/01_reproducible_workflow.ipynb |
| 测试肽已有 Cys3/Cys14 评分 | results/2mi1_reduced_capped_ph75_20260924/analysis/results.csv 和 notebooks/02_reduced_peptide_case.ipynb |
| 结构文件与关联 | 同案例 start_*/ 中的 PDB；CSV structure_file 是仓库相对路径 |
| 第三方版本、参数、调用 | docs/THIRD_PARTY.md、logs/reduced_md_20260924、protocol.json、requirements.txt |
| 数据处理、划分与泄漏 | data/README.md、provenance；明确训练片段与测试肽重叠，不当作独立验证 |
| 打分、排序、不确定性 | 原始分数与过滤结果保存；每帧内排序，不解释为效率；短采样无收敛声明 |
| 新序列插入及后续验证 | design.py / app.py 支持本地 ESMFold 自动结构预测与评分；KTTKS 演示已保存；尚未完成 GFN2-xTB 能垒或湿实验验证 |
| 默认结果字段 | 案例清单含 ID、赛道占位、序列/修饰、位点、分数、结构、模型和运行版本、备注 |
| 本地复核包 | tools/package_submission.py 输出 ZIP 与逐文件 SHA-256，无上传动作 |

旧版模型在还原态案例中未稳定复现实验排序，不将该案例包装为成功验证。赛事反馈说暂无专用模板，本项目采用默认 CSV 字段；实际赛道仍需填写。原始实验材料、数据许可和全面同源/预训练重叠核验仍待补齐，不能因完成整理而宣称这些材料齐备。


当前正式训练为 logs/retrain_sst_20260924，SST 已作为一条实验排序观察纳入训练；结果见 results/sst_ranking_training_20260924 与 [SST 重训说明](SST_TRAINING.md)。旧版还原态和 PROPKA 报告仅作历史对照。
