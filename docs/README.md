# 文档索引

项目介绍、统一安装命令和评审演示见 [根目录 README](../README.md)。所有 Python 依赖仅在 [requirements.txt](../requirements.txt) 维护，覆盖基础运行、Notebook、ESMFold、OpenMM 和 PROPKA；本目录中的命令默认在仓库根目录执行。

| 阅读目的 | 文档 |
| --- | --- |
| 陌生 Windows 电脑双击安装与启动 | [一键启动](ONE_CLICK_START.md) |
| 网页使用、本地结构预测与评估包导出 | [自动设计](AUTOMATIC_DESIGN.md) |
| 插入位置、手工结构输入与筛选规则 | [插入设计](DESIGN_WORKFLOW.md) |
| 湿实验反馈字段、审核及训练划分 | [反馈接口](FEEDBACK.md) |
| 外部能垒请求、反馈与独立回归模型 | [能垒接口](BARRIER_INTERFACE.md) |
| 当前正式模型的排序监督与复现 | [SST 重训说明](SST_TRAINING.md) |
| 旧版模型的还原态模拟与 PROPKA 对照 | [还原态案例](REDUCED_PEPTIDE_CASE.md) |
| KTTKS 文献依据与候选准备记录 | [KTTKS 案例](KTTKS_CASE_STUDY.md) |
| 正式目录、历史归档与打包范围 | [项目结构](PROJECT_LAYOUT.md) |
| 赛事交付材料及尚待补充的信息 | [提交核对表](SUBMISSION_CHECKLIST.md) |
| 第三方用途、版本与来源记录 | [第三方说明](THIRD_PARTY.md)、[文件来源哈希](file_origins.json) |

数据与结果的就地说明分别见 [数据说明](../data/README.md)、[结果说明](../results/README.md)和 [Model Card](../models/MODEL_CARD.md)。`logs/` 下的执行快照及历史日志保留当时内容，复现以当前入口和本目录指南为准。

`docs/` 是纳入提交包的正式指南；`tmp/docs/` 保存赛事原件、反馈和本地参考资料，`tmp/archive/` 保存历史版本与过程文件。它们不属于提交内容，正式运行不依赖其中的文件。ESMFold 权重首次运行时单独下载；重做历史 MD 模拟另需配置系统 OpenCL 驱动。
