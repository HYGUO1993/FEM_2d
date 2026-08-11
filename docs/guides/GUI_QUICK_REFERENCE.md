# FEM_2d GUI 快速参考

## 快速启动

- Windows 双击：
  - `scripts/run_gui.bat`
  - `scripts/run_gui.ps1`
  - 或根目录 `quick_start_gui.py`
- 命令行：
  - PowerShell: `.\scripts\run_gui.ps1`
  - CMD: `cmd /c scripts\run_gui.bat`
  - Python: `python quick_start_gui.py`
  - Linux/macOS: `python3 python/fem_gui.py`

## 环境要求

- Python 3.8+
- 依赖包：`PyQt6 matplotlib numpy`
- 安装命令：
  - `pip install PyQt6 matplotlib numpy`
  - 或 `conda install -c conda-forge PyQt6 matplotlib numpy`

## 常用操作流程

1. 打开模型（`File > Open` 或示例文件）
2. 调整显示参数（变形倍率、反力倍率）
3. 运行求解（`Analysis > Solve` 或 `F5`）
4. 浏览结果（缩放/平移/显示开关）
5. 导出图像（`Analysis > Export Plot...`）

## 常用快捷键

- `Ctrl+O`：打开模型
- `Ctrl+S`：保存模型
- `F5`：运行求解
- `F6`：切换变形显示
- `F7`：切换轴力着色
- `F8`：切换反力显示

## 常见问题

- 缺失包：`python -m pip install PyQt6 matplotlib numpy`
- 点击求解无可执行文件：先运行 `bash scripts/build_and_test.sh`（或 PowerShell 版本）
- DLL 加载失败：见 `docs/guides/WINDOWS_DLL_ERROR.md`

## 相关文档

- `README.md`
- `docs/guides/GUI_USAGE.md`
- `docs/guides/USAGE_GUIDE.md`
- `docs/development/DEVELOPMENT.md`
- `RELEASE_NOTES.md`
