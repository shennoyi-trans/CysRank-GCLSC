# Cys 插入候选结构报告

枚举 6 个插入位置，评分 6 个，输出 3 条不同序列。

优先选择通过过滤的候选；不足时按内部评分补足并提示风险。各组内按分数降序，序列去重。模型分数不是成功概率。

结构检查仅含重原子完整性、主链键长与严重非局部重叠，不能证明折叠正确或反应可行。历史训练特征口径不完整，存在分布偏移风险。

## 候选 1

序列：`KTTCKS`
在原序列第 3 位之后插入 C；新位点编号 4。
模型分数：0.962518；0.5 阈值分类：True。
Cys 侧链 SASA：89.651 Å²；SG 周围 5 Å 非氢原子数：8。
结构来源：ESMFold v1 / Hugging Face Transformers。结构文件：structures/peptide_insert_C_after_3.pdb。

## 候选 2

序列：`KTTKCS`
在原序列第 4 位之后插入 C；新位点编号 5。
模型分数：0.945512；0.5 阈值分类：True。
Cys 侧链 SASA：89.198 Å²；SG 周围 5 Å 非氢原子数：8。
结构来源：ESMFold v1 / Hugging Face Transformers。结构文件：structures/peptide_insert_C_after_4.pdb。

## 候选 3

序列：`KCTTKS`
在原序列第 1 位之后插入 C；新位点编号 2。
模型分数：0.942074；0.5 阈值分类：True。
Cys 侧链 SASA：97.057 Å²；SG 周围 5 Å 非氢原子数：8。
结构来源：ESMFold v1 / Hugging Face Transformers。结构文件：structures/peptide_insert_C_after_1.pdb。
