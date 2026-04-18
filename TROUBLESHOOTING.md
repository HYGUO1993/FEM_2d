# 📋 FEM_2d 故障排除工具清单

## 🔍 诊断工具

### `python diagnose.py`
**功能**: 完整的环境诊断
- ✅ 检查 Python 版本和路径
- ✅ 验证所有依赖包（PySide6, matplotlib, numpy）
- ✅ 检查 FEM 求解器编译状态
- ✅ **精确识别 Windows DLL 加载错误**
- ✅ 提供针对性的解决方案

**何时使用**: 任何问题出现时，首先运行此工具

**示例输出**:
```
[!!] WINDOWS DLL LOADING ERROR DETECTED
Solutions (try in order):
  1. Install Visual C++ Redistributable
  2. Reinstall PySide6 after installing VC++
  3. Try conda install
```

---

## 🔧 自动修复工具

### `修复DLL错误.bat` / `fix_dll_error.bat`
**功能**: 一键修复 Windows DLL 加载错误 ⭐

**做什么**:
1. 检测 Visual C++ 运行时是否已安装
2. 如果缺失，自动下载并安装
3. 清除 pip 缓存
4. 重新安装 PySide6
5. 验证修复是否成功

**何时使用**: 看到 "DLL load failed" 或 "找不到指定的程序" 错误

**使用方法**: 
- 双击运行（推荐）
- 或在 PowerShell 运行: `./fix_dll_error.bat`

**系统要求**:
- Windows 系统
- 网络连接（自动下载 Visual C++ 运行时）
- 如果自动安装失败，需要管理员权限手动安装

---

### `python quick_start_gui.py`
**功能**: 智能依赖检查和启动

**做什么**:
1. 检查所有依赖包
2. 如果缺失，尝试自动安装
3. 提供多种安装方法备选
4. 最终启动 GUI

**何时使用**: 依赖缺失或需要启动 GUI

**使用方法**:
```bash
python quick_start_gui.py
```

---

## 🚀 启动工具

### `启动GUI.bat` / `start_gui.bat`
**功能**: 启动 FEM_2d 图形界面

**何时使用**: 诊断和修复完成后，启动 GUI

**使用方法**: 双击或 `启动GUI.bat`

---

## 📚 文档和参考

### `FIX_GUI_ERRORS.md`
**内容**: 快速错误修复指南
- DLL 加载错误的 3 种修复方法
- 缺失包的安装方法
- 无法启动时的命令行替代方案

**何时查看**: 看到特定错误信息时

---

### `WINDOWS_DLL_ERROR.md`
**内容**: Windows DLL 问题详细解决方案
- 问题原因详解
- 6 种修复方法
- Visual C++ 运行时安装指南
- conda 替代方案
- 故障排除步骤

**何时查看**: DLL 错误的详细信息或自动修复失败

---

### `QUICKSTART_Windows.md`
**内容**: Windows 专用快速启动指南
- 最简单的启动方式
- 分步操作指南
- 快捷方式表
- 常见问题解答
- 命令行使用说明

**何时查看**: 第一次使用或需要详细指导

---

## 🔄 典型故障排除流程

```
1. 遇到问题
   ↓
2. 运行 python diagnose.py
   ↓
3. 根据诊断结果选择修复工具
   ├─→ DLL 错误 → 运行 修复DLL错误.bat
   ├─→ 缺失包 → 运行 python quick_start_gui.py
   └─→ 其他 → 查看 FIX_GUI_ERRORS.md
   ↓
4. 验证修复: python diagnose.py
   ↓
5. 启动 GUI: 启动GUI.bat 或 python quick_start_gui.py
```

---

## 📞 如果仍然有问题

1. **获取诊断信息**
   ```bash
   python diagnose.py > diagnostics.txt
   ```

2. **查看相关文档**
   - `FIX_GUI_ERRORS.md` - 快速修复指南
   - `WINDOWS_DLL_ERROR.md` - DLL 问题详解
   - `QUICKSTART_Windows.md` - Windows 快速指南

3. **尝试命令行方式**
   ```bash
   build\bin\Release\fem_run.exe --input test05.txt
   ```

---

## 🎯 工具快速参考

| 问题 | 诊断工具 | 修复工具 | 文档 |
|------|---------|---------|------|
| GUI 无法启动 | `diagnose.py` | `修复DLL错误.bat` | `FIX_GUI_ERRORS.md` |
| DLL 加载失败 | `diagnose.py` | `修复DLL错误.bat` | `WINDOWS_DLL_ERROR.md` |
| 缺失依赖包 | `diagnose.py` | `quick_start_gui.py` | `INSTALL_DEPENDENCIES.md` |
| 首次使用 | - | - | `QUICKSTART_Windows.md` |

---

**最后更新**: 2026年4月18日

✨ **提示**: 大多数问题可以通过 `诊断 → 修复 → 启动` 三步解决！
