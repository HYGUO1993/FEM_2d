# FEM_2d 快速启动指南 (Windows)

## 🚀 最简单的方式（推荐）

### **如果出现 "DLL load failed" 错误：**

**双击运行修复工具：**
```
修复DLL错误.bat
```

或（英文版）：
```
fix_dll_error.bat
```

这个工具会自动：
1. ✅ 检测 Visual C++ 运行时是否已安装
2. ✅ 如果缺失，自动下载并安装
3. ✅ 重新安装 PySide6
4. ✅ 验证修复是否成功

**然后启动 GUI：**
```
启动GUI.bat
```

---

## 📋 分步指南

### 步骤 1: 检查问题
```bash
python diagnose.py
```

看到什么错误？

- **"DLL load failed"** → 使用 `修复DLL错误.bat`
- **"No module named PySide6"** → 运行 `python quick_start_gui.py`
- **其他错误** → 查看 `FIX_GUI_ERRORS.md`

### 步骤 2: 运行对应的修复

#### 如果是 DLL 错误（最常见）
```bash
修复DLL错误.bat
```

#### 如果是缺失包
```bash
python quick_start_gui.py
```

#### 或使用 Conda（最可靠）
```bash
conda install -c conda-forge PySide6 matplotlib numpy
```

### 步骤 3: 启动 GUI
```bash
启动GUI.bat
```

或
```bash
python quick_start_gui.py
```

---

## 🎯 快捷方式

| 目的 | 命令 | 说明 |
|------|------|------|
| 修复 DLL 错误 | `修复DLL错误.bat` | 自动安装 Visual C++ + 重装 PySide6 |
| 启动 GUI | `启动GUI.bat` | 直接打开图形界面 |
| 诊断问题 | `python diagnose.py` | 查看详细的环境信息 |
| 快速安装依赖 | `python quick_start_gui.py` | 自动检查并安装依赖 |

---

## ❓ 常见问题

### Q: 运行 `修复DLL错误.bat` 后仍然有问题？

A: 尝试这些步骤：
1. 重启计算机（安装 VC++ 后需要重启）
2. 运行 `python diagnose.py` 检查当前状态
3. 查看 `WINDOWS_DLL_ERROR.md` 获取其他解决方案
4. 尝试使用 Conda：`conda install -c conda-forge PySide6`

### Q: 安装 Visual C++ 需要多久？

A: 通常 2-5 分钟，取决于网络速度和系统。

### Q: 可以使用命令行替代 GUI 吗？

A: 可以！查看 `USAGE_GUIDE.md` 了解命令行用法。

---

## 📁 主要文件说明

| 文件 | 用途 |
|------|------|
| `修复DLL错误.bat` / `fix_dll_error.bat` | 一键修复 Windows DLL 错误 |
| `fix_dll_error.ps1` | PowerShell 修复脚本（自动执行） |
| `启动GUI.bat` / `start_gui.bat` | 启动图形界面 |
| `quick_start_gui.py` | Python 依赖检查和启动工具 |
| `diagnose.py` | 环境诊断工具 |
| `FIX_GUI_ERRORS.md` | 快速错误修复指南 |
| `WINDOWS_DLL_ERROR.md` | Windows DLL 问题详细解决方案 |

---

## 🔧 命令行使用

如果 GUI 无法运行，可以使用命令行：

```bash
# 运行求解器
build\bin\Release\fem_run.exe --input test05.txt --output Results.dat

# 生成可视化
python scripts/visualize_results.py --results Results.dat --out plot.png

# 查看结果
start plot.png
```

---

**最后更新**: 2026年4月18日

💡 **提示**: 如果一切正常，应该在运行 `启动GUI.bat` 后看到图形界面。
