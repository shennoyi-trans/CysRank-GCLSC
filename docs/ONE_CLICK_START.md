# Windows 一键启动

将整个项目解压到有写权限的本地目录，双击根目录 `start.cmd`。无需先激活 Conda、输入 pip 命令或手动打开浏览器。支持 Windows x64（Intel/AMD）；当前锁定依赖对应 CPython 3.13 普通 64 位版本。首次安装需要能访问 Python.org、PyPI、PyTorch 下载站和 Hugging Face，并有足够磁盘空间和内存。ESMFold 权重约 8.44 GB，Python、依赖和安装缓存还需额外空间；没有足够 GPU 显存时使用 CPU，计算可能较慢。

## 自动执行的步骤

1. 检查项目环境、PATH、Python Launcher、Python 注册表及 Conda 环境列表，选择兼容且可复用依赖最多的解释器，优先复用体积较大的 Torch。
2. 没有兼容解释器时，从 [Python 官方发布页](https://www.python.org/downloads/release/python-3135/)下载已验证的 3.13.5 x64 安装包，检查 Python Software Foundation 数字签名，按[官方静默安装接口](https://docs.python.org/3.13/using/windows.html#installing-without-ui)安装到项目 `.runtime/python313`，不修改 PATH。
3. 已有项目虚拟环境直接复用。外部解释器在项目 `.runtime/venv-313` 建立可继承已有包的虚拟环境；缺失或版本不符的包安装到项目环境，避免修改外部环境。复制后失效的旧环境会保留，自动另建新环境。此策略复用的是本机已有环境，不能把虚拟环境跨电脑直接搬运。
4. 使用唯一的 `requirements.txt` 补齐完整依赖，版本满足时跳过安装；保留 pip 缓存。执行依赖一致性检查，并实际导入 Torch、NumPy、Biopython、Transformers、Accelerate 与项目评分模块。
5. 优先复用项目 `models/folding-cache`，其次复用标准 Hugging Face 缓存（含 `HF_HOME`/`HF_HUB_CACHE` 设置）。没有完整缓存时，按项目固定模型修订下载配置及权重。遵循 [Hugging Face 缓存下载接口](https://huggingface.co/docs/huggingface_hub/guides/download)，不重复下载已经完整缓存的文件。准备完成后后端使用本地模型路径，预测时无需再下载。
6. 启动后端，默认绑定 `127.0.0.1:8765`；端口被占用时自动选空闲端口。在服务绑定成功后打开默认浏览器。启动窗口会显示实际地址；关闭窗口或按 Ctrl+C 可停止服务。

同一目录的启动器持有排他锁，防止重复双击导致并发安装。已有服务运行时，请使用已打开的页面；需要重新启动时先关闭原启动窗口。模型和依赖准备失败不会打开一个貌似可用的页面；错误会保留在窗口和 `logs/launcher/*.log` 中，修复网络或文件问题后重新双击即可。

## 排查与命令行

在项目目录运行以下命令可仅检查环境和缓存，不安装、不下载、不启动后端：

```powershell
.\start.cmd -CheckOnly
```

需要使用指定端口或不打开浏览器时：

```powershell
.\start.cmd -Port 8877 -NoBrowser
```

启动器不是离线依赖包：陌生电脑首次运行需要联网，或预先准备匹配的完整环境和模型缓存。系统管理策略禁用脚本执行、网络不可达、磁盘不足、缺少完整项目文件等情况会明确报错。完整 ESMFold 推理仍需足够系统内存，启动检查不代替每条序列的模型推理验证。

## 交付

`start.cmd`、`tools/start.ps1`、`tools/launch.py` 及本说明均纳入项目打包范围。`.runtime`、虚拟环境、启动日志及模型缓存不纳入提交包；新电脑会重新探测、复用该电脑现有内容并安装缺失部分。仅复制启动脚本不能运行，必须连同源码、`requirements.txt`、`web`、`docs` 及项目随附权重完整交付。

## 本次验证（2026-10-03）

- Windows 本机通过 `start.cmd` 选择可复用依赖最多的解释器，创建 `.runtime/venv-313`，复用已有 Torch 等包，补装缺失包，未重新下载 ESMFold。
- 再次启动跳过安装，默认浏览器自动请求页面与配置接口，均返回 HTTP 200。
- 启动器 9 项测试、网页后端 3 项测试和已有结构流程 2 项测试通过，覆盖缺失依赖、间接依赖修复、缓存命中/缺失、环境复用、保留失效环境、端口冲突及浏览器启动。
- 新后端真实执行 `CC → CCC`：三个插入位置完成评分，输出一条去重候选；外部评估 ZIP 经 HTTP 下载并通过完整性校验。任务目录为 `results/web_jobs/41bd928abbe640dbb316a144fcd0481a/`。该极短肽用于工程验证，结果为补足候选，不是实验效果证明。
- 打包文件清单已确认包含启动入口、辅助脚本与说明。没有在另一台完全未安装 Python 的电脑实测；自动下载安装 Python 的分支采用官方安装接口和签名校验。
