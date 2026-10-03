# 第三方软件与调用记录

| 工具 | 本项目用途与版本 | 来源与许可说明 |
| --- | --- | --- |
| ProteinMPNN | 冻结 CA-only 编码器，v_48_020.pt；原始代码及权重 SHA-256 见 file_origins.json | https://github.com/dauparas/ProteinMPNN；MIT 文本保留在 licenses/ProteinMPNN-MIT.txt；原始下载时间/commit 未提供 |
| OpenMM | 8.6.1；2026-09-24 本地还原态 MD，Amber ff14SB/TIP3P，OpenCL | https://github.com/openmm/openmm；各组件主要采用 MIT/LGPL，OpenCL 后端为 LGPL，见官方许可链接 |
| Biopython | 1.88；PDB 读取及 Shrake–Rupley SASA | https://biopython.org；随包许可适用 |
| PROPKA | 3.5.1；固定 λ=0.5 的 pKa 后处理对照 | https://github.com/jensengroup/propka；按上游许可使用；通过 requirements.txt 安装 |
| PDBFixer | 1.12.0；仅用于未完成的前期修复排查，正式封端模拟未调用 | https://github.com/openmm/pdbfixer；MIT；不作为正式运行依赖 |
| ESMFold v1 | 新增 Cys 候选结构预测；facebook/esmfold_v1，revision 75a3841ee059df2bf4d56688166c8fb459ddd97a | https://github.com/facebookresearch/esm；MIT；权重从 https://huggingface.co/facebook/esmfold_v1 下载，仅缓存于本地 |
| Transformers / Accelerate | 4.57.6 / 1.12.0；本地 ESMFold 推理与低内存加载 | https://github.com/huggingface/transformers / https://github.com/huggingface/accelerate；Apache-2.0；运行库通过 requirements.txt 安装，不复制入提交包 |

OpenMM 官方安装说明：https://docs.openmm.org/latest/userguide/application/01_getting_started.html

OpenMM 许可说明：https://docs.openmm.org/latest/userguide/library/01_introduction.html ；许可原文：https://github.com/openmm/openmm/blob/master/docs-source/licenses/Licenses.txt 。正式包不复制这些第三方运行库，通过锁定版本的依赖文件安装。

所有正式运行依赖统一在 requirements.txt 中声明，包括 OpenMM、Transformers、Accelerate 和 PROPKA；仅重算已有结构的特征/评分不会调用 OpenMM。软件版本、设备、随机种子、步长、力场、溶剂和质子化假设都记录在案例 protocol.json 和 logs/reduced_md_20260924；各模拟 system.xml 保存实际参数化系统，原执行脚本按运行哈希留存。

本项目没有调用外部 API 或商业结构预测服务。原始序列、结构及实验标签的来源许可按 data/README.md 披露，不由软件许可证替代。数据和项目自身代码的发布授权须由对应权利人确认。
