# 自动结构预测与多肽输入窗口

本流程完成：输入陌生多肽 → 枚举单 Cys 插入 → 本地 ESMFold 结构预测 → 现有 CysRank 评分 → Top3 → 外部 GFN2-xTB 交接包。原有正式分类头不重训。

## 启动

在仓库根目录安装唯一的完整依赖清单，然后启动服务：

```powershell
python -m pip install -r requirements.txt
python app.py
```

浏览器打开 http://127.0.0.1:8765 。保留目标多肽输入框，可反复输入新序列。输入自动去除空白并转大写；只接受标准氨基酸，至少 2 个残基。网页默认上限 100 个残基，可通过 `--max-length` 调整；上限是任务规模限制，不是精度保证。

已安装依赖且模型缓存完整时，可直接执行 `python app.py --local-files-only`。环境配置见 [项目快速开始](../README.md#快速开始与评审演示)。

首次使用下载 facebook/esmfold_v1 权重（约 8.44 GB）至 models/folding-cache。序列推理在本地执行，不上传到结构预测服务。后续可用 `--local-files-only` 离线运行。`--fold-model` 可指定本地 Hugging Face ESMFold 模型目录。

`--device auto` 仅在可用 GPU 显存至少 12 GiB 时选择 CUDA，否则使用 CPU；阈值是保守配置，不是显存保证。ESM 语言模型部分使用 float16，结构模块保留 float32。分块默认 32，可通过 `--chunk-size` 调整。CPU 仍需要足够系统内存，长序列和大量插入候选可能耗时较长。

网页同一服务一次处理一个预测任务，避免多个大模型同时占用内存。刷新页面可以重新连接本次任务；结果存于 results/web_jobs/<任务编号>/design，服务重启后磁盘结果仍保留。

## 命令行

```powershell
python design.py run --sequence KTTKS --record-id example --output results/example_auto --local-files-only
```

输出目录必须是新路径。输出包含 preparation（全部插入映射、去重序列和预测结构）、results（评分及报告）、external_evaluation.zip（外部评估包）。已有 `prepare` / `score` 命令继续支持人工提供结构。

## 排名规则

1. 枚举 N+1 个插入操作；每条不同序列预测一次结构。
2. 对每个插入操作的新增 Cys 提取特征并评分，保留已有过滤和结构诊断结果。
3. 先取通过过滤的候选，按原始模型分数降序排列。
4. 不足三条时，从未通过过滤但有有限数值分数的候选中按分数补足。
5. 最终按完整序列去重。相同序列优先保留通过过滤的解释。同组同分时以 candidate_id 稳定排序。

补足项标记 `selection_status=fallback`，并在网页、CSV、JSONL 和报告中提示“根据已有了解，该序列的成功可能性较低”，附具体原因。该提示是基于规则的定性判断，不是概率估计。几何异常也明确提示需修复。缺结构、序列不匹配、缺重原子或非有限分数不能用虚构分数补足。

若不同序列本身不足三条（例如 CC 的所有插入操作都得到 CCC），输出实际数量，并说明原因。原 `eligible` 仍表示是否通过过滤，不因进入 Top3 而被改成 true。

## 外部评估包

包中包含 top3.csv、top3.jsonl、全候选评分、重原子 PDB、结构预测来源、模型与输入哈希、报告及能垒请求。external_barrier_requests.jsonl 使用包内相对结构引用，保持 sample 和 sample_sha256 不变。

GFN2-xTB 在外部执行。交接的肽 PDB 不是完整反应计算体系：计算方需要定义试剂、反应步骤、补氢/端基、质子化、电荷、溶剂、参考态、路径与过渡态核验。未计算的能垒为 null，不将结构优化能或反应能作为能垒。详细字段见 BARRIER_INTERFACE.md。

## 验证和范围

2026-09-25：40 项自动测试通过；AG 的三条插入候选通过真实本地 ESMFold 完成预测、CysRank 评分和 ZIP 导出。此极短肽的预测结构触发几何警告，已按规则提示。

此外，从网页输入 KTTKS，已完成全部六个插入位置的真实预测、评分与下载 ZIP 校验。Top3 为 KTTCKS（0.9625）、KTTKCS（0.9455）、KCTTKS（0.9421），均通过当前过滤。结果保存于 results/kttks_esmfold_20260925。这些是工程运行结果，不是反应质量或模型泛化验证；没有执行 GFN2-xTB。

第一版每条候选使用一个未松弛的预测构象；尚未加入 MD、多构象聚合或反应建模。结构置信度不作为实验成功概率。项目内不存在 local-preview 目录，网页为 Python 服务直接提供的静态 HTML，无 npm 构建依赖。

官方模型：https://huggingface.co/facebook/esmfold_v1

实现接口：https://huggingface.co/docs/transformers/v4.57.6/model_doc/esm

ESM 上游：https://github.com/facebookresearch/esm （MIT）；Transformers 及其 OpenFold 移植部分遵循各自上游许可。模型缓存和网页运行任务不纳入提交包；打包脚本包含输入界面、自动预测模块和安装说明。
