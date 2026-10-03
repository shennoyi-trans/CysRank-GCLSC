# 2MI1 还原态已有位点复核

本案例回答反馈意见（3）中的具体问题：对同一条测试肽的已有 Cys3 和 Cys14 评分，不插入新残基。实验排序由提供方说明为 Cys3 环化效率明显高于 Cys14。本案例使用旧版 `part2-20260924-v1` 模型，其计算排序未稳定复现该观察；当前正式模型已纳入 SST 排序监督，见 [重训说明](SST_TRAINING.md)。

## 实验说明与建模边界

2026-09-24 用户确认：反应前已还原，氧化态不能反应；N 端乙酰化 Ac、C 端一级酰胺 NH2、Cys 巯基无保护、三氟乙酸盐、pH 7.5、300 K。反应缓冲体系及 TFA 浓度未知，计算采用 0.15 M NaCl 显式水近似，未显式模拟 TFA。

原输入 data/structures/test.pdb 为含 10 个构象的 2MI1。它的 SG–SG 连接及约 2 Å 距离对应氧化态，但同时含两个 HG，连接与氢信息不一致。构建时删除旧氢及 SG–SG 键，增加 ACE/NHE 端基并重新补氢，两个 Cys 显式使用中性巯基；不等于恒 pH 模拟或位点 pKa 计算。

## 实际流程与结果

OpenMM 8.6.1，Amber ff14SB/TIP3P，NVIDIA RTX 4050 的 OpenCL 后端，NVT，300 K，2 fs 步长。每个起点最小化后运行 20 ps 平衡和 100 ps 生产采样，每 10 ps 留一帧，共 100 个生产快照。整轮实测约 347 秒，不含安装和评分。10 个起点使用不同固定随机种子，未按分数选择构象。

| 结果 | 数值 |
| --- | --- |
| Cys3 平均原始分数 | 0.891827 |
| Cys14 平均原始分数 | 0.985284 |
| Cys3 分数更高的生产快照 | 7/100 |
| 按每个起点平均后 Cys3 更高 | 0/10 |
| 生产快照 SG–SG 距离 | 3.45–11.14 Å |
| 当前几何检查/过滤失败 | 0/100 |

结果不是环化效率；100 帧有时间相关性，不是 100 次独立实验。100 ps 采样不能证明收敛。模型训练含 AGCKNFFW 和 KTFTSC，与测试肽相关，因此不作为独立验证。SASA/邻近原子数包含 Ac/NH2 的重原子，但模型序列编码仅使用 14 个标准残基。与前一轮仅修改二硫键标记的试验相比，端基和构象同时改变，不能单独归因于还原。

## 查看与复算

- [案例报告](../results/2mi1_reduced_capped_ph75_20260924/analysis/REPORT.md)：全流程与局限。
- [标准化逐位点清单](../results/2mi1_reduced_capped_ph75_20260924/analysis/results.csv)：200 行生产快照位点结果，含序列、修饰、位点、分数、模型版本/哈希及对应 PDB。每帧内排序；不将两百行称为新候选分子。
- [全部结构比较](../results/2mi1_reduced_capped_ph75_20260924/analysis/comparison.csv)：另含最小化、平衡结构。
- [可执行案例 Notebook](../notebooks/02_reduced_peptide_case.ipynb)：显示实际结果并重新提取特征、评分核验。
- results/2mi1_reduced_capped_ph75_20260924/start_*/：完整带端基结构，单位 Å，链 A；原肽残基 1–14，ACE=0、NHE=15。
- logs/reduced_md_20260924/start_*/：系统/积分器 XML、溶剂化初始结构、最终状态/检查点、能量轨迹；executed_scripts/ 保存原实际执行脚本，仅用于溯源。

安装统一依赖后复算随附结构的特征和分数（此步骤不调用 OpenMM、不训练模型）：

```powershell
python -m pip install -r requirements.txt
python tools/analyze_2mi1_reduced_md.py --output results/2mi1_reanalysis
```

重新生成模拟结构（同样使用 requirements.txt 中的 OpenMM；另需 OpenCL 驱动，输出目录须不存在）：

```powershell
python -m pip install -r requirements.txt
python tools/rebuild_2mi1_reduced.py --output results/2mi1_new_run --log-output logs/2mi1_new_run --ph 7.5 --temperature 300 --opencl-platform 0 --device 0
python tools/analyze_2mi1_reduced_md.py --input results/2mi1_new_run
```

OpenCL 平台和设备编号需按机器调整；本机平台 0 为 NVIDIA。不同驱动/GPU 上轨迹不保证逐坐标一致，应核验配置、拓扑、能量稳定性和统计结果。导出的二进制检查点也具有平台依赖性，XML 和初始结构一并保留。

早期准备失败、pH 7 中断运行、仅改标记的控制试验完整移到 tmp/archive/reduced_case_process_20260924；不参与正式案例汇总。正式推理、模型与训练数据均未修改。


PROPKA 固定 λ=0.5 的探索性对照已完成：Cys3 更高的生产快照由 7/100 变为 6/100，未改善实验排序。包含端基类型适配及其局限，见 [完整报告](../results/2mi1_propka_lambda05_20260924/REPORT.md)。
