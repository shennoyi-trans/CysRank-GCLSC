# CysRank SST 排序监督模型

# SST 排序监督重训

版本 `sst-ranking-20260924-v1`；50 epochs、seed=20260916，AdamW lr=0.01、weight_decay=0.01。固定训练设置，未按排序结果追加训练或选权重。

保留 19 条分类记录、22 个监督位点（8 正、14 负）；新增一条用户提供的实验排序 Cys3 > Cys14。100 个还原封端生产构象只作该观察的输入增广，不把 Cys14 标成阴性。原始氧化 PDB 本身不提供反应标签。实验原始测量未提供。

损失：`(sum BCE + mean_conformers softplus(-(z3-z14)))/23`。分类每位点权重为 1；整条 SST 排序观察权重为 1，平均到 100 帧，每帧占 0.01。不是 100 个额外实验。PROPKA 修正不参与训练，λ 未修改。

| 指标 | 结果 |
| --- | --- |
| 分类训练回代 | 22/22 |
| SST 训练构象 Cys3 > Cys14 | 99/100（旧版 7/100） |
| Cys3 平均分 | 0.983623 |
| Cys14 平均分 | 0.769741 |
| 每轨迹平均后 Cys3 更高 | 10/10 |
| 全 SST 来源组留出后的排序 | 26/100 构象 |

分组诊断中 SST、AGCKNFFW、KTFTSC 均属于 old3_positive_05；留出该组时所有排序构象和相关片段同时移出训练。该诊断说明移除 SST 监督后的表现，不是训练后模型的新样本验证；100 帧仍非独立样本。其余诊断缺少独立阴性，不能声称总体泛化准确率。

本次改进说明拟合到了已知实验排序，不证明新肽筛选有效。后续应冻结模型，使用新的来源组、预先确定的候选和对照做统一能垒/实验验证。端基未进入序列编码、短时 MD 未收敛、历史负例完整蛋白环境等局限仍存在。

## 复现

```powershell
python train.py --output logs/sst_reproduction
python run.py --output logs/sst_workflow
```

configs/train.json 已默认启用 data/ranking/sst_ordering.json；删除配置中的 ranking 键可复现纯分类训练（另用新输出目录）。正式权重为 models/site_predictor.pt。历史 MD/PROPKA 对照仍绑定 logs/retrain_20260924/final_model/site_predictor.pt，不会被新版权重悄然替换。旧版完整快照在 tmp/archive/before_sst_training_20260924。

逐构象新旧分数见 results/sst_ranking_training_20260924/comparison.csv；训练、分组消融、原始数据快照、哈希和损失见 logs/retrain_sst_20260924。

模型 SHA-256：`58f1c0ec37d005d971078b6d6d148f6b7fe693f7d5c4951db9d67c699629b1d4`
