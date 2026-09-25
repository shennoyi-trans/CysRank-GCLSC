# 还原态案例运行证据

protocol.json、runs.json、snapshots.jsonl 为原运行清单的逐字节副本。executed_scripts 为当时真实运行的源码快照，可与 protocol/analysis 中的 script_sha256 核对；不要在此目录直接执行这些历史文件，正式复现入口位于 tools。

start_01 至 start_10 保存系统与积分器 XML、显式溶剂初始结构、最终状态/二进制检查点及 energy.csv。每条轨迹 60,000 步（20 ps 平衡 + 100 ps 生产），不是完整逐步轨迹；10 ps 间隔的肽结构保存在 results/2mi1_reduced_capped_ph75_20260924。

preliminary_control 保存此前仅修改二硫键标记的对照汇总，彼时端基及坐标未重建，不纳入最终模拟统计。失败/中断尝试保留在 tmp/archive/reduced_case_process_20260924。
