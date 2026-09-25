# 外部 GFN2-xTB 评估交接

top3.csv / top3.jsonl 为模型排序；warning 标记补足候选及原因。分数不是实验成功概率。

structures/ 为预测肽的重原子结构。external_barrier_requests.jsonl 的 structure_reference 相对本文件；sample 保留原始模型输入用于哈希追溯。不要修改 sample 或 sample_sha256。

计算方需补氢、设定端基/质子化、电荷、溶剂、试剂及明确反应步骤，构建反应物/产物及路径，核验过渡态后报告 ΔE‡ 或 ΔG‡。本项目未执行 xTB，未填写能垒值。

详细字段见 BARRIER_INTERFACE.md。几何检查失败的结构需先修复；失败计算记录不能填为零。
