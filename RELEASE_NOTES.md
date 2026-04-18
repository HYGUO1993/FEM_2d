# FEM_2d Release v1.0

## 概述

FEM_2d 是一个用于教学和工程实践的二维杆系/桁架/刚架有限元分析程序。本版本包含完整的C++核心计算引擎和Python可视化脚本。

## 主要功能

- ✅ 二维桁架（Truss）和刚架（Frame）结构分析
- ✅ 支持多种材料和截面类型
- ✅ 节点集中力载荷
- ✅ 位移、内力和支座反力计算
- ✅ Python 可视化（变形图、内力图、支座反力）

## 快速开始

### 1. 构建程序

**Linux/macOS:**
```bash
./scripts/build_and_test.sh
```

**Windows (PowerShell):**
```powershell
./scripts/build_and_test.ps1
```

或手动构建：
```bash
mkdir build && cd build
cmake ..
cmake --build .
```

### 2. 运行示例

```bash
# 悬臂梁示例
./build/bin/fem_run --input test05.txt --output build/Results.dat

# 简支梁示例
./build/bin/fem_run --input test_beam.txt --output build/Results_beam.dat
```

### 3. 可视化结果

```bash
# 安装 Python 依赖
pip install matplotlib numpy

# 生成可视化图像
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --out plot.png
```

## 文件说明

### 核心文件
- `barsystem.cpp/h` - FEM 核心计算引擎
- `test05.txt` - 悬臂梁示例
- `test_beam.txt` - 简支梁示例

### 脚本文件
- `scripts/visualize_results.py` - 结果可视化脚本
- `scripts/build_and_test.sh` - Linux/macOS构建脚本
- `scripts/build_and_test.ps1` - Windows构建脚本

### GUI 应用 (Beta)
- `python/fem_gui.py` - 图形化界面入口（需要 PySide6）
- `python/gui/` - GUI 模块

## 输入文件格式

第一行为总控数据：
```
节点数 约束节点数 单元数 材料类型数 截面类型数 载荷数
```

示例（test05.txt）：
```
2 1 1 1 1 1
2 0.0 0.0
2 1.0 0.0
0 -1 -1 -1
2 0 1 0
210000000000 0.3 0.0
0.01 8.333e-6 0.1
1 1 -1000 -1 1 0 0 0
```

详细格式说明请参见 README.md

## 命令行参数

- `--input <file>` - 输入文件路径（默认: test05.txt）
- `--output <file>` - 输出文件路径（默认: Results.dat）
- `--quiet` - 静默模式
- `--stiff <file>` - 输出单元刚度矩阵
- `--no-stiff` - 不输出单元刚度矩阵

## 可视化选项

- `--results <file>` - 结果文件路径
- `--scale <float>` - 变形倍率
- `--out <file>` - 输出图片路径
- `--hide-reactions` - 隐藏支座反力
- `--hide-forces` - 隐藏内力着色
- `--reac-scale <float>` - 反力箭头比例

## 系统要求

**编译环境:**
- CMake 3.10+
- C++11 编译器（GCC 7+, Clang 5+, MSVC 2017+）

**可视化:**
- Python 3.7+
- matplotlib
- numpy

**GUI (可选):**
- PySide6
- PyVista

## 测试

运行单元测试：
```bash
./build/bin/unit_tests
```

或使用 CTest：
```bash
cd build
ctest
```

## 已知限制

- 当前仅支持节点集中力载荷
- 不支持分布载荷（均布、三角形等）
- 不支持温度载荷
- GUI 功能仍在开发中

## 开发文档

详细的开发文档请参见：
- `README.md` - 完整使用说明
- `DEVELOPMENT.md` - 开发者文档
- `.github/copilot-instructions.md` - 代码规范

## 版本历史

### v1.0 (2026-04-18)
- 初始发布版本
- 完整的 FEM 核心功能
- Python 可视化脚本
- 基础 GUI 框架
- 单元测试

## 许可证

请参见项目仓库的 LICENSE 文件

## 联系方式

- GitHub: https://github.com/HYGUO1993/FEM_2d
- Issues: https://github.com/HYGUO1993/FEM_2d/issues

---

**构建信息：**
- 分支: claude/improve-visualization-for-engineers
- 编译日期: 2026-04-18
- 核心代码: C++11
- 测试状态: ✅ 通过
