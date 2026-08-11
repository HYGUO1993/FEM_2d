#!/usr/bin/env python3
"""
FEM_2d GUI Quick Reference Guide

This file documents all available methods to launch and use the FEM_2d GUI application.
For detailed information, see GUI_USAGE.md
"""

# ============================================================
# QUICK START - 快速开始
# ============================================================

"""
Method 1: Double-click (Windows) - 双击启动 ⭐ 推荐
─────────────────────────────────────────────────────
1. Open File Explorer 打开文件浏览器
2. Navigate to: scripts/ 导航到 scripts 文件夹
3. Double-click one of:
   - run_gui.bat        (Batch script / 批处理脚本)
   - run_gui.ps1        (PowerShell script / PowerShell 脚本)
   - Or in project root: quick_start_gui.py
"""

# Method 2: Command Line - 命令行启动
"""
Windows PowerShell:
  .\scripts\run_gui.ps1
  
Windows CMD:
  cmd /c scripts\run_gui.bat
  
Python direct:
  python quick_start_gui.py
  
macOS/Linux:
  python3 python/fem_gui.py
  ./scripts/run_gui.ps1
"""

# Method 3: IDE Integration - 集成到 IDE
"""
VS Code:
  1. Open project in VS Code
  2. Run task: Tasks > Run Task > Run GUI
     (Requires .vscode/tasks.json configuration)

PyCharm:
  1. Right-click quick_start_gui.py
  2. Select "Run 'quick_start_gui'"
  
CLion:
  1. Create Run Configuration
  2. Script: quick_start_gui.py
  3. Python interpreter: (select your environment)
"""

# ============================================================
# ENVIRONMENT SETUP - 环境配置
# ============================================================

"""
1. Python Requirements 最低要求
   - Python 3.8 or higher
   - Conda (optional but recommended)

2. Install Dependencies 安装依赖
   Manual installation:
   ──────────────────
   pip install PySide6 matplotlib numpy
   
   Using conda:
   ──────────
   conda install -c conda-forge PySide6 matplotlib numpy
   
   Automatic (via launch scripts):
   ─────────────────────────────
   Launch scripts will automatically install missing packages

3. Build Release Solver (Optional but Recommended)
   ───────────────────────────────────────────────
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --config Release
   
   This is needed for the "Solve" (F5) button to work
"""

# ============================================================
# FEATURES - 功能说明
# ============================================================

"""
GUI Features 图形界面功能:
───────────────────────────
✓ Interactive visualization      交互式可视化
✓ Real-time parameter control    实时参数控制
✓ Open/save model files          打开/保存模型
✓ One-click analysis             一键分析
✓ Export results to images       导出结果图像
✓ Model tree view                模型树视图
✓ Properties panel               属性面板
✓ Run logs                        运行日志
✓ Deformation scaling            变形缩放
✓ Reaction force visualization   反力可视化
✓ Axial force coloring           轴力着色
"""

# ============================================================
# BASIC WORKFLOW - 基本工作流程
# ============================================================

"""
Typical usage flow 典型工作流程:
─────────────────────────────

1. OPEN MODEL 打开模型
   File > Open  (Ctrl+O)
   - Load your .txt model file
   
   OR
   
   File > Open Example (test05)
   - Use built-in example

2. ADJUST PARAMETERS (可选)
   - On the right panel:
   - Deformation scale
   - Reaction scale
   - Show/hide display options

3. RUN SOLVER 运行求解器
   Analysis > Solve  (F5)
   - Watch the logs for progress
   - Visualization updates automatically

4. EXPLORE RESULTS 查看结果
   - Pan/zoom with mouse
   - Use toolbar buttons for navigation
   - Adjust display options on right panel

5. EXPORT RESULTS (可选)
   Analysis > Export Plot...
   - Save as PNG, JPEG, or SVG
"""

# ============================================================
# KEYBOARD SHORTCUTS - 键盘快捷键
# ============================================================

"""
File Operations:
  Ctrl+N  - New model
  Ctrl+O  - Open model
  Ctrl+S  - Save model
  Ctrl+Q  - Exit

Analysis:
  F5      - Solve
  
View Control:
  F6      - Toggle deformed shape
  F7      - Toggle force coloring
  F8      - Toggle reactions
"""

# ============================================================
# TROUBLESHOOTING - 故障排除
# ============================================================

"""
Issue: "Python not found"
─────────────────────────
Solution:
  1. Ensure Python is installed
  2. Add Python to PATH
  3. Use conda environment

Issue: "Module not found" (PySide6, matplotlib, numpy)
──────────────────────────────────────────────────────
Solution:
  pip install PySide6 matplotlib numpy
  
Issue: "Solver not found" when clicking Solve
────────────────────────────────────────────
Solution:
  Build Release version:
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build --config Release

Issue: Visualization not updating
────────────────────────────────
Solution:
  1. Click Redraw button
  2. Check logs for error messages
  3. Verify model file format

See GUI_USAGE.md for more detailed troubleshooting
"""

# ============================================================
# DOCUMENTATION - 文档
# ============================================================

"""
Main Documents 主要文档:
─────────────
  README.md              - Project overview
  GUI_USAGE.md           - Complete GUI guide ⭐ READ THIS
  USAGE_GUIDE.md         - Command-line tool guide
  DEVELOPMENT.md         - Developer documentation
  RELEASE_NOTES.md       - Version history

API Reference 代码参考:
─────────────
  python/fem_gui.py              - Application entry point
  python/gui/main_window.py       - Main GUI window
  python/gui/__init__.py          - Module initialization
"""

# ============================================================
# FILE LOCATIONS - 文件位置
# ============================================================

"""
Project Structure 项目结构:
──────────────────────────

FEM_2d/
├── python/                    # Python GUI code
│   ├── fem_gui.py            # Application entry
│   ├── test_api.py           # API tests
│   └── gui/
│       ├── __init__.py       # Package init
│       └── main_window.py    # GUI window class
│
├── scripts/                   # Helper scripts
│   ├── run_gui.bat           # Windows batch launcher
│   ├── run_gui.ps1           # PowerShell launcher
│   ├── run_visualize.bat     # CLI visualization
│   └── build_and_test.sh/ps1 # Build scripts
│
├── build/                     # Build output
│   ├── bin/Release/
│   │   └── fem_run.exe       # Solver executable
│   ├── Results_gui.dat       # Last analysis result
│   └── plot_gui.png          # Last generated plot
│
├── GUI_USAGE.md              # ⭐ Complete GUI guide
├── quick_start_gui.py        # Quick launcher
├── test05.txt                # Example model
└── test_beam.txt             # Example model
"""

# ============================================================
# VERSION INFORMATION - 版本信息
# ============================================================

"""
FEM_2d GUI
Version: 1.0.0
Last Updated: April 18, 2026

For latest version and releases:
https://github.com/HYGUO1993/FEM_2d/releases
"""

if __name__ == "__main__":
    print(__doc__)
