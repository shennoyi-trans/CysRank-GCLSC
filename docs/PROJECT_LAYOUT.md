# 正式交付目录与过程归档

2026-09-24 依据附件4与反馈意见（3）整理；2026-10-03 统一依赖清单并恢复根目录 docs 正式文档。模型、原始 PDB、训练数据和正式模拟快照保留原内容并校验哈希。

| 目录 | 正式内容 |
| --- | --- |
| configs、data、src、models | 训练配置、数据与来源、核心算法、正式权重及模型说明 |
| notebooks | 01 训练复现；02 还原态测试肽案例及重算核验 |
| results | 标准训练回代示例、KTTKS 插入准备材料与本地 ESMFold 评分、pH 7.5 还原态结构及旧模型评分 |
| logs/retrain_20260924 | 旧版分类模型训练、诊断及重载证据 |
| logs/retrain_sst_20260924 | 当前正式分类与排序联合训练、诊断及重载证据 |
| logs/verification_20260924 | 本轮模型测试与 notebook 验证 |
| logs/reduced_md_20260924 | 模拟参数、种子、系统、能量、检查点及原执行脚本 |
| logs/organization_20260924 | 本次重算、Notebook、打包与迁移核验 |
| docs、licenses、tools、tests | 正式指南、第三方说明、CLI 工具及回归测试 |

tmp/docs 保留赛事要求原件、反馈和历史说明；正式使用指南以根目录 docs 为准。tmp/agents 是本机临时依赖、辅助脚本和验证副本。tmp/archive/reduced_case_process_20260924 保存此次失败/中断尝试及旧日志，MOVES.json 逐文件记录移动和 SHA-256。

正式工具不读取 tmp 中的运行时；基础运行、Notebook、ESMFold、MD 和 PROPKA 均按唯一的 requirements.txt 配置环境，重做 MD 另需系统 OpenCL 驱动。打包排除 tmp、.git、编辑器缓存、本机环境和重复 notebook 运行目录；logs 中保留的是必要复现证据。执行脚本的历史快照仅用于核对原哈希，复现使用 tools 中的当前入口。


当前正式训练为 logs/retrain_sst_20260924，SST 已作为一条实验排序观察纳入训练；结果见 results/sst_ranking_training_20260924 与 [SST 重训说明](SST_TRAINING.md)。旧版还原态和 PROPKA 报告仅作历史对照。
