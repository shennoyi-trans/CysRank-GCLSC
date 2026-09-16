# 项目整理与恢复记录

整理日期：2026-09-16。以原 result_from_ZQY(8) 为本次交付成果，并从其他目录提取其必需的训练源码、ProteinMPNN 权重、MIT 许可及母蛋白分组数据。最终模型和原始数据保持字节不变；23 个复制文件的源路径和 SHA-256 见 file_origins.json。

## 目录映射

| 原位置 | 正式项目位置 |
| --- | --- |
| ZQY(8) 正负训练文件 | data/raw |
| ZQY(8) 验证文件 | data/external/validation.jsonl |
| ZQY(8) 删除前正例备份 | data/provenance/positive_before_exclusion.jsonl |
| ZQY(3) 第二部分正例 | data/provenance/parent_sequences.jsonl，仅用于分组 |
| ZQY(8) 测试.pdb | data/structures/test.pdb |
| ZQY(8) final_model/site_predictor.pt | models/site_predictor.pt |
| ZQY(8) prepared/training.jsonl | data/processed/training.jsonl |
| ZQY(8) 原配置、训练报告、日志和诊断 | logs/original |
| ca_reproducibility 的必要代码 | src/core.py、src/pipeline.py、src/protein_mpnn_utils.py |
| ca_reproducibility 的官方 CA-only 权重 | models/pretrained/v_48_020.pt |
| 附件4及来源方要求 | docs/submission_requirements.docx、source_requirements.docx |

训练实现仅调整模块组织、相对导入、版本标识及路径记录，移除不再使用的旧版准备和硬过滤函数；没有修改训练特征、优化器或最终模型。新增根目录训练/预测/筛选/完整流程入口，CSV 导出、路径解析和文档。源代码保留在归档中以供对照。

## 归档与恢复

所有旧成果版本、旧代码目录、第三方完整仓库、旧文档、Obsidian 配置、渲染缓存、本机 .venv-gpu 及 Python 缓存统一移到 `tmp/archive_before_submission_20260916/`，未删除。精确移动清单位于该目录的 moves.json。根目录 `.git` 为版本管理元数据，保留原位，不参与提交包。

file_origins.json 中 source 是整理前的仓库相对位置；现在可以在上述归档目录下找到对应原件。恢复时按 moves.json 的 relative_source 从归档移回，先检查目标是否存在，避免覆盖正式项目。无需恢复旧目录即可运行当前项目。

tmp/portability_check 保存移植验证副本，tmp 中其他脚本为本次整理辅助文件，均不应进入评审包。评审只需正式项目文件和新建的 Python 环境，不依赖归档虚拟环境。

## 验证

12 项回归测试覆盖数据清洗、监督掩码、反馈隔离及闭环、标签不进入模型特征和结构平移不变性。23 份源文件复制哈希通过；原成果 19 个位点与当前模型重新推理分数差为 0。复制至新目录并从外部工作目录启动预测，CSV 字节一致，训练流程也通过。完整 50 轮训练、优化和 Top 3 导出已运行，结果保存在 logs/verification。

torch.cross 的弃用警告来自随附第三方原始代码；本次保留其实现和许可证，没有修改第三方计算逻辑。原始日志中的历史绝对路径和机器信息作为证据保留；运行代码与配置不依赖它们。
