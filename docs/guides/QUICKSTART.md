# FEM_2d GUI 启动指南

> **最快 3 秒启动 GUI！**

---

## 选择您的启动方式

### ⭐ 方式 1：双击启动（Windows 用户推荐）

最简单、最方便的方式！

1. **打开项目文件夹**
   - 找到 `FEM_2d` 项目根目录

2. **双击以下文件之一**
   - `启动GUI.bat` （中文用户）
   - `start_gui.bat` （英文用户）

3. **等待 GUI 窗口弹出** ✓
   - 第一次启动会自动安装依赖（需要网络）
   - 后续启动速度很快

> ✅ 一键启动，无需输入命令

---

### 🔌 方式 2：VS Code 集成（编程用户）

在 VS Code 中直接启动，无需打开其他窗口

**快捷键启动：**
1. 按 `Ctrl+Shift+P`
2. 输入 `Launch FEM_2d GUI`
3. 按 `Enter`

**或通过菜单：**
- 终端 → 运行任务 → 🚀 Launch FEM_2d GUI

> ✅ 在 VS Code 中无缝集成，编码分析两不误

---

### 💻 方式 3：命令行启动

适合高级用户或脚本集成

```bash
# 最推荐：跨平台（Windows/macOS/Linux）
python quick_start_gui.py

# 或直接调用
python python/fem_gui.py

# Windows 批处理
.\start_gui.bat
.\scripts\run_gui.bat
```

> ✅ 灵活控制，支持参数和脚本集成

---

## 首次启动

首次启动时，脚本会：

1. ✓ 检查 Python 环境
2. ✓ 自动检查缺失的依赖包（PySide6、matplotlib、numpy）
3. ✓ 自动安装缺失的包（如需要）
4. ✓ 检查 FEM 求解器是否可用
5. ✓ 启动 GUI 应用

**需要：** 网络连接（用于安装依赖）

---

## 遇到问题？

### 问题：双击脚本后闪退

**原因：** Python 未安装或环境问题

**解决：**
1. 确保 Python 已安装（使用 `python --version` 验证）
2. 如使用 Anaconda，先激活环境：`conda activate base`
3. 尝试改用方式 3（命令行启动）

### 问题：缺少依赖包

**原因：** 自动安装失败或无网络

**解决：**
```bash
pip install PySide6 matplotlib numpy
```

### 问题："Solve" 按钮无法使用

**原因：** FEM Release 求解器未编译

**解决：**
```bash
# 在项目根目录执行
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

### 更多问题？

查看完整指南：[GUI_USAGE.md](GUI_USAGE.md)

---

## 基本操作（3 步分析）

1. **打开模型** 
   - `File` → `Open` 或 `Open Example`

2. **求解** 
   - `Analysis` → `Solve`（或按 `F5`）

3. **查看结果** 
   - 在中央画布区交互查看
   - 右侧面板调整显示参数

---

## 键盘快捷键速查

| 快捷键 | 功能 |
|--------|------|
| `Ctrl+O` | 打开模型 |
| `F5` | 求解 |
| `Ctrl+Q` | 退出 |

更多快捷键见 [GUI_USAGE.md](GUI_USAGE.md)

---

## 文件说明

### 启动脚本
| 文件 | 说明 |
|------|------|
| `启动GUI.bat` | ⭐ Windows 双击启动（中文） |
| `start_gui.bat` | Windows 双击启动（英文） |
| `quick_start_gui.py` | Python 跨平台启动 |
| `scripts/run_gui.bat` | 完整版启动脚本 |
| `scripts/run_gui.ps1` | PowerShell 启动脚本 |

### 文档
| 文件 | 说明 |
|------|------|
| `GUI_USAGE.md` | 📖 完整 GUI 使用手册 |
| `GUI_QUICK_REFERENCE.md` | 快速参考卡 |
| `README.md` | 项目说明 |

---

## 下一步

✅ 已启动 GUI？

→ 查看 [GUI_USAGE.md](GUI_USAGE.md) 学习完整功能

---

**版本：** 1.0.0  
**更新：** 2026年4月18日  
**许可证：** MIT
