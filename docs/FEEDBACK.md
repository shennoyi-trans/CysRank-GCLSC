# 反馈接口

每行 JSON 需要 observation_id、evidence_type（wet_experiment / simulation）、split（train / holdout）、parent_group、assay_id、非空 conditions 对象、非空 provenance、sample（完整 seq/xyz/四项生化特征）、position_1based、outcome（0/1/null）、reviewed（布尔值）。可选 cyclization_efficiency_fraction 为 [0,1] 比例，不接受百分数，也不从效率自动推断二元结局。

调用方式见 README。输入文件需来源方自行提供；本次没有虚构湿实验反馈。归档仅校验格式，不代表自动审核实验真实性。训练只接收已审核的 train 湿实验记录，模拟、留出和未审核记录只归档。重复位点和母蛋白分组冲突需人工处理，不能反复提高样本权重。
