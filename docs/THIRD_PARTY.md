# 第三方说明

ProteinMPNN 原仓库标识 https://github.com/dauparas/ProteinMPNN，本地随附 MIT 许可（Copyright 2022 Justas Dauparas），全文保留于 licenses/ProteinMPNN-MIT.txt。未提供下载时间、Git commit 或下载日志，因此以随附文件 SHA-256 标识本地快照，不编造上游提交号。

src/protein_mpnn_utils.py 原样复制自现有运行版本；来源和哈希见 file_origins.json。CA-only 预训练权重 v_48_020.pt 的 SHA-256 为 f28f40170e21858c5ff31ef50b6e63414ff76dc331b19f85aa8586a12031744a。本地权重随 MIT 仓库保留，原下载来源及是否有专属附加条款需权利人进一步核对。

调用方式为 Python 本地实例化：CA-only、128 隐藏维、3 层 encoder/decoder、48 邻居、augment_eps=0、dropout=0。训练记录日期 2026-09-16，输入为坐标、序列及结构特征；详见 logs/original/config.json。不使用外部 API、商业平台、提示词或在线结构生成服务。

Python、PyTorch、NumPy 及传递依赖版本见 requirements.txt 和 logs/verification/environment.txt，分别遵守各自许可证。项目自身代码及实验数据的发布授权应由对应权利人确认。
