# KTTKS：有文献依据的陌生线性多肽设计案例

检索与本地核对日期：2026-09-23。以下保留当日的候选准备与文献核对记录。2026-09-25 已完成六个插入位置的本地 ESMFold 结构预测与评分，见 [自动设计结果](../results/kttks_esmfold_20260925/results/report.md)；尚无能垒或改造后的实验结果。

## 选择依据

母序列为 **KTTKS（Lys–Thr–Thr–Lys–Ser）**，本案例采用未脂化、全 L 型、游离 N 端和羧酸 C 端的线性形式 H-KTTKS-OH。这是设计输入的明确约定；实际采购、建模和实验时须核对端基、盐型及质子化状态。它属于有护肤研究背景的肽，不将其称为已获批透皮药物。

- Katayama 等，1993，*A pentapeptide from type I procollagen promotes extracellular matrix production*：细胞研究将 KTTKS 确定为促进胶原和纤连蛋白生成的最小活性序列。[原始论文摘要，PMID 8486721](https://pubmed.ncbi.nlm.nih.gov/8486721/)。这是细胞层面的依据，不等于改造后活性或人体抗皱效果已获证实。
- Choi 等，2014，*Dermal Stability and In Vitro Skin Permeation of Collagen Pentapeptides (KTTKS and palmitoyl-KTTKS)*：在其离体无毛小鼠皮肤实验条件下，KTTKS 未在各皮层检出，两种肽均未在接收液检出；棕榈酰化版本在各皮层有检出且更稳定。[原始论文，DOI 10.4062/biomolther.2014.053](https://pmc.ncbi.nlm.nih.gov/articles/PMC4131521/)。该结果支持研究局部皮肤递送瓶颈，不应外推为任何配方、任何物种下均无法吸收。

选择未修饰 KTTKS 可直接使用现有标准氨基酸输入流程；pal-KTTKS 的棕榈酰基不在当前结构/模型支持范围，不能只输入 KTTKS 就声称评价了 pal-KTTKS。

## 本地陌生性核对

以母序列及六条插入序列为查询，对 data 下现有所有 JSONL 的 seq/sequence 做精确相等及双向完整子串检查，匹配记录均为 0：

| 文件 | 记录数 | 匹配记录数 |
|---|---:|---:|
| raw/训练正例.jsonl | 5 | 0 |
| raw/训练集负例.jsonl | 11 | 0 |
| processed/training.jsonl | 16 | 0 |
| provenance/parent_sequences.jsonl | 16 | 0 |
| provenance/positive_before_exclusion.jsonl | 6 | 0 |
| external/validation.jsonl | 250 | 0 |
| examples/input.jsonl | 16 | 0 |

文件间存在重复，不应把记录数相加作为独立样本数。这只证明当前列明数据未发现上述序列重叠，不证明全球首次设计、无同源性或 ProteinMPNN 预训练从未见过相关结构。

## 已生成的候选

运行命令：`python design.py prepare --sequence KTTKS --record-id kttks_literature_case --output results/kttks_design_requests`。

| 新序列 | 新 Cys 位置（1-based） | 2026-09-23 准备状态 |
|---|---:|---|
| CKTTKS | 1 | 保留枚举记录；游离 N 端不满足当前相邻主链酰胺规则，评分时排除 |
| KCTTKS | 2 | 待结构与反应评估 |
| KTCTKS | 3 | 待结构与反应评估 |
| KTTCKS | 4 | 待结构与反应评估 |
| KTTKCS | 5 | 待结构与反应评估 |
| KTTKSC | 6 | 待结构与反应评估 |

KTTKSC 保留了连续的 KTTKS 序列，可作为优先讨论的实验假设；添加 Cys 同时改变了原 Ser 的端基环境，仍可能改变活性，不能据此认定最佳。其余内部插入会打断原活性序列，更需检验功能保留。该准备目录不包含分数或 Top3；后续评分保存在独立的 `results/kttks_esmfold_20260925/`。

输出目录含 candidates.jsonl、candidate_sequences.fasta、structure_manifest.template.json 和 preparation.json。评分需要各候选的配套全原子结构；原始 KTTKS 的结构不能替代插入后的结构。短肽的构象多样性还需结构准备人员处理，单个构象排名可能不稳健。

## 本案例要检验的假设

项目所用 [Sun 团队 2025 年反应](https://doi.org/10.1021/jacs.5c08837)是 Cys 侧链与主链参与的局部成环；不是将整条肽首尾闭合。具体结构与条件参见该论文及其补充材料；本地参考副本保存在 tmp/docs/references/Sun_2025_JACS/，不纳入提交包。当前没有证据说明该反应必然提高 KTTKS 的皮肤递送，CA 修饰本身也可能影响性质。

建议把目标写为：**在保留功能活性的前提下，探索 Cys 插入及 CA 局部环化对 KTTKS 稳定性和皮肤内递送的影响。** 对护肤应用应区分表皮/真皮到达量与穿过全层皮肤进入接收液的量。

建议的最小对照为原始 KTTKS、对应插入但未环化的六肽、经结构确认的 CA 环化产物；有条件可另加入文献中的 pal-KTTKS 作为递送参照。比较应控制剂量、载体、时间和检测方法，分别测定产物身份、稳定性、皮层分布及功能活性。若要进一步区分 CA 加成与成环的贡献，还需要合适的化学对照。

能垒接口只服务于反应可行性筛选，不输出皮肤通透性或活性标签。此案例仍是待验证设计，不可作为已成功改善透皮的赛事结果。
