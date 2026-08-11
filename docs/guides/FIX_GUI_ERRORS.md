# 🔴 GUI 启动失败？快速修复！

## 诊断问题

首先，运行诊断工具来精确定位问题：

```bash
python diagnose.py
```

这会显示环境状态和具体错误信息。

---

## 根据错误类型修复

### 情况 1: "DLL load failed" 或 "找不到指定的程序"

这是 **Windows 系统库缺失**，不是 Python 包问题。

**快速修复（3 步）：**

1. 下载并安装 Visual C++ 运行时：
   - 访问：https://support.microsoft.com/en-us/help/2977003
   - 下载：Visual C++ Redistributable（x64 版本）
   - 双击安装，重启计算机

2. 重新安装 PySide6：
   ```bash
   python -m pip cache purge
   python -m pip install --force-reinstall PySide6
   ```

3. 尝试启动 GUI：
   ```bash
   python quick_start_gui.py
   ```

**或使用 conda（更可靠）：**
```bash
conda install -c conda-forge PySide6 matplotlib numpy
```

👉 详细步骤见：[WINDOWS_DLL_ERROR.md](WINDOWS_DLL_ERROR.md)

---

### 情况 2: "No module named PySide6"

这是 Python 包**未安装**。

**修复：**
```bash
python -m pip install --user PySide6 matplotlib numpy
```

或：
```bash
conda install -c conda-forge PySide6 matplotlib numpy
```

---

### 情况 3: 其他错误

运行诊断并查看错误信息，然后查看相关文档：
- 📖 [GUI_USAGE.md](GUI_USAGE.md) - GUI 完整使用手册
- 📖 [INSTALL_DEPENDENCIES.md](INSTALL_DEPENDENCIES.md) - 依赖问题详解
- 📖 [WINDOWS_DLL_ERROR.md](WINDOWS_DLL_ERROR.md) - Windows 特定问题

---

## 无法启动 GUI 时的替代方案

如果 GUI 始终无法工作，可以使用命令行：

```bash
# 运行求解器
./build/bin/Release/fem_run.exe --input test05.txt --output build/Results.dat

# 生成可视化
python scripts/visualize_results.py --results build/Results.dat --out build/plot.png

# 查看结果
start build/plot.png
```

---

## 还有问题？

1. 再次运行诊断：`python diagnose.py`
2. 根据输出选择相应的解决方案
3. 查看详细的文档指南
4. 或使用命令行工具替代方案

---

**最后更新**：2026年4月18日
