# Windows DLL 加载错误修复

## 问题症状

```
Error: DLL load failed while importing QtWidgets: 找不到指定的程序
```

或

```
Error: The procedure entry point could not be located in the dynamic link library
```

## 原因

PySide6（PyQt）在 Windows 上需要以下系统库：
- Visual C++ 运行时库
- OpenGL 支持
- Windows Media Foundation

这些是系统级依赖，仅安装 Python 包是不够的。

---

## 解决方案

### 方案 1：安装 Visual C++ 可再发行组件包（最有效！）

从 Microsoft 官网下载并安装：

**Visual Studio 2015-2022 可再发行组件包**  
https://support.microsoft.com/en-us/help/2977003

下载 "Visual C++ 可再发行程序包" → 选择 x64 版本（如果 Python 是 64 位）

安装后重启计算机。

---

### 方案 2：安装完整的 Visual Studio 构建工具（推荐开发者）

```
https://visualstudio.microsoft.com/visual-cpp-build-tools/
```

选择"Desktop development with C++"组件。

---

### 方案 3：重新安装 PySide6（有时有效）

```bash
# 清理缓存
python -m pip cache purge

# 强制重新安装
python -m pip install --force-reinstall --no-cache-dir PySide6

# 或卸载后重新安装
python -m pip uninstall -y PySide6
python -m pip install PySide6
```

---

### 方案 4：使用 PyQt5 替代（绕过方案）

如果 PySide6 仍然无法工作，可以改用 PyQt5：

1. **修改 `python/gui/main_window.py`**：
   - 将 `from PySide6.QtWidgets import ...` 改为 `from PyQt5.QtWidgets import ...`
   - 将 `from PySide6.QtCore import ...` 改为 `from PyQt5.QtCore import ...`
   - 将 `from PySide6.QtGui import ...` 改为 `from PyQt5.QtGui import ...`

2. **修改 `python/fem_gui.py`**：
   - 同样改为 `from PyQt5...`

3. **安装 PyQt5**：
   ```bash
   python -m pip install PyQt5
   ```

4. **修改 matplotlib 后端**（在 `python/gui/main_window.py` 开头）：
   ```python
   import matplotlib
   matplotlib.use('Qt5Agg')  # 改为 Qt5 后端
   ```

---

### 方案 5：在虚拟环境中使用 Conda（最可靠）

如果问题持续存在，使用 conda 创建干净的环境可能更有效：

```bash
# 创建新的 conda 环境
conda create -n fem2d python=3.10 -y

# 激活环境
conda activate fem2d

# 安装依赖（conda 处理系统库更好）
conda install -c conda-forge PySide6 matplotlib numpy -y

# 运行 GUI
python python/fem_gui.py
```

---

### 方案 6：命令行工具替代（完全绕过 GUI）

如果 GUI 无法工作，可以使用命令行：

```bash
# 运行求解器
./build/bin/Release/fem_run.exe --input test05.txt --output build/Results.dat

# 生成可视化（无需 GUI）
python scripts/visualize_results.py --results build/Results.dat --out build/plot.png

# 打开生成的图像
start build/plot.png
```

---

## 诊断步骤

1. **检查诊断报告**：
   ```bash
   python diagnose.py
   ```
   查看"IMPORT TEST"部分的错误信息

2. **查看 DLL 依赖**（高级用户）：
   ```bash
   # 使用 dependency walker 检查 PySide6 的 DLL 依赖
   # 下载：depends.exe 或使用 Python 工具
   python -c "import PySide6; print(PySide6.__path__)"
   ```

3. **重启计算机**：
   某些系统库安装后需要重启才能生效

---

## 测试修复

安装上述库后，测试是否有效：

```bash
# 运行诊断
python diagnose.py

# 尝试启动 GUI
python quick_start_gui.py
```

如果诊断中"IMPORT TEST"的 PySide6 行现在显示 ✓，说明问题已解决。

---

## 如果仍然无法工作

1. 查看"IMPORT TEST"的具体错误信息
2. 搜索错误信息到 Stack Overflow 或 PySide 官方文档
3. 考虑使用虚拟机或 WSL (Windows Subsystem for Linux)
4. 或使用命令行替代方案（见方案 6）

---

## 参考资源

- PySide6 官方文档：https://doc.qt.io/qtforpython/
- PyQt5 官网：https://www.riverbankcomputing.com/software/pyqt/
- Visual C++ 运行时：https://support.microsoft.com/en-us/help/2977003
- Qt 官方：https://www.qt.io/

---

**最后更新**：2026年4月18日
