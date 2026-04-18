# FEM_2d v1.0 使用指南

## 快速开始 - 3步即可使用

### 步骤1：获取代码
```bash
git clone https://github.com/HYGUO1993/FEM_2d.git
cd FEM_2d
git checkout claude/improve-visualization-for-engineers
```

### 步骤2：编译程序
```bash
# Linux/macOS
./scripts/build_and_test.sh

# 或者手动编译
mkdir build && cd build
cmake ..
make
```

### 步骤3：运行示例
```bash
# 运行悬臂梁分析
./build/bin/fem_run --input test05.txt --output build/Results.dat

# 查看结果
head -50 build/Results.dat
```

---

## 功能特性

### ✅ 已实现功能

1. **结构类型**
   - 二维桁架（Truss）
   - 二维刚架（Frame）

2. **分析功能**
   - 线性静力分析
   - 节点位移计算
   - 单元内力计算
   - 支座反力计算

3. **载荷类型**
   - 节点集中力（X、Y、弯矩方向）

4. **求解器**
   - Skyline 存储格式
   - LDLT 分解求解

5. **可视化**
   - 结构变形图
   - 内力着色显示
   - 支座反力箭头
   - 可调变形倍率

---

## 使用示例

### 示例1：悬臂梁

**输入文件 (test05.txt):**
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

**说明：**
- 2个节点，1个约束节点，1个单元
- 节点0固定，节点1自由端
- 节点1施加竖向-1000N力

**运行：**
```bash
./build/bin/fem_run --input test05.txt --output Results_cantilever.dat
```

### 示例2：简支梁

**运行：**
```bash
./build/bin/fem_run --input test_beam.txt --output Results_beam.dat
```

---

## 可视化

### 安装依赖
```bash
pip install matplotlib numpy
```

### 生成图片
```bash
# 基本用法
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --out plot.png

# 隐藏反力箭头
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --hide-reactions --out plot.png

# 隐藏内力着色
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --hide-forces --out plot.png

# 自定义反力箭头比例
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --reac-scale 0.1 --out plot.png
```

---

## 输入文件格式详解

### 第1行：总控数据
```
节点总数 约束节点数 单元数 材料类型数 截面类型数 载荷数
```

### 节点数据（NodeType X Y）
```
2 0.0 0.0    # 节点类型2（Frame），坐标(0, 0)
2 1.0 0.0    # 节点类型2（Frame），坐标(1, 0)
```
- NodeType: 1=桁架节点，2=刚架节点

### 约束数据（Node DX DY DR）
```
0 -1 -1 -1   # 节点0，XYR方向全约束
```
- -1 = 约束该自由度
- 0 = 自由

### 单元数据（Type Node1 Node2 Section Material）
```
2 0 1 0 0    # Frame单元，连接节点0和1，使用截面0和材料0
```
- Type: 1=桁架，2=刚架

### 材料数据（E μ α）
```
210000000000 0.3 0.0
```
- E: 弹性模量（Pa）
- μ: 泊松比
- α: 热膨胀系数

### 截面数据（A Iz H）
```
0.01 8.333e-6 0.1
```
- A: 面积（m²）
- Iz: 惯性矩（m⁴）
- H: 截面高度（m）

### 载荷数据（Type Dir Value Elem Node Pos T0 T1）
```
1 1 -1000 -1 1 0 0 0
```
- Type: 1=节点集中力
- Dir: 0=X, 1=Y, 2=R（转角）
- Value: 载荷值
- Node: 节点编号
- Elem: 单元编号（-1表示不适用）

---

## 输出文件说明

### 节点位移 (Node Displacements)
```
Node        Ux        Uy        Rz
   0  0.000000  0.000000  0.000000
   1  0.001905 -0.001905  0.002857
```

### 单元端力 (Element End Forces)
```
Elem    Fx_i    Fy_i    Mz_i    Fx_j    Fy_j    Mz_j
   0       0    1000     500       0   -1000       0
```
- i端：单元起点
- j端：单元终点

### 支座反力 (Support Reactions)
```
Node        Rx        Ry        Rz
   0         0     -1000      -500
```

---

## 图形界面（GUI）- Beta版

### 安装GUI依赖
```bash
pip install PySide6 pyvista
```

### 启动GUI
```bash
python python/fem_gui.py
```

### GUI功能
- 📁 文件菜单：新建、打开、保存
- ✏️ 编辑菜单：撤销、重做
- 🔬 分析菜单：求解
- 👁️ 视图菜单：显示变形、显示内力
- 🌳 模型树：查看节点、单元、载荷等
- 📊 属性面板：查看对象属性

**注意：** GUI功能仍在开发中，部分功能尚未完全集成。

---

## 命令行参数完整列表

### fem_run 参数
```
--input <file>     输入文件路径（默认: test05.txt）
--output <file>    输出文件路径（默认: Results.dat）
--quiet            静默模式，减少输出
--stiff <file>     输出单元刚度矩阵到文件
--no-stiff         不输出单元刚度矩阵
```

### visualize_results.py 参数
```
--results <file>   结果文件路径（默认: build/Results.dat）
--scale <float>    变形倍率（默认: 1.0）
--out <file>       输出图片路径（默认: build/plot.png）
--hide-reactions   隐藏支座反力箭头
--hide-forces      隐藏单元内力着色
--reac-scale <f>   反力箭头比例（默认: 自动）
```

---

## 常见问题

### Q: 编译失败？
A: 确保安装了CMake 3.10+和C++11编译器：
```bash
# Ubuntu/Debian
sudo apt-get install cmake g++

# macOS
brew install cmake
```

### Q: Python可视化报错？
A: 安装依赖：
```bash
pip install matplotlib numpy
```

### Q: 结果位移为0？
A: 检查是否使用了正确的单元类型：
- 梁/刚架问题必须使用Frame单元（类型2）
- 确保截面惯性矩Iz > 0

### Q: 位移过大或NaN？
A: 检查约束是否充分：
- 平面结构至少需要3个约束自由度
- 确保结构不会发生刚体位移

### Q: GUI无法启动？
A: 安装PySide6：
```bash
pip install PySide6
```

---

## 文件清单

```
FEM_2d/
├── barsystem.cpp          # 核心计算引擎
├── barsystem.h            # 头文件
├── CMakeLists.txt         # CMake配置
├── test05.txt             # 悬臂梁示例
├── test_beam.txt          # 简支梁示例
├── README.md              # 项目说明
├── RELEASE_NOTES.md       # 发布说明
├── USAGE_GUIDE.md         # 本文件
├── DEVELOPMENT.md         # 开发文档
├── scripts/
│   ├── visualize_results.py    # 可视化脚本
│   ├── build_and_test.sh       # Linux/macOS构建
│   └── build_and_test.ps1      # Windows构建
├── python/
│   ├── fem_gui.py         # GUI入口
│   └── gui/               # GUI模块
│       ├── __init__.py
│       └── main_window.py
├── tests/
│   └── test_basic.cpp     # 单元测试
└── build/                 # 构建目录
    └── bin/
        ├── fem_run        # 主程序
        └── unit_tests     # 测试程序
```

---

## 下一步学习

1. **阅读完整文档：** README.md
2. **查看开发文档：** DEVELOPMENT.md
3. **运行所有示例：** test05.txt 和 test_beam.txt
4. **创建自己的模型：** 参考输入文件格式
5. **尝试GUI：** python python/fem_gui.py

---

## 获取帮助

- **GitHub**: https://github.com/HYGUO1993/FEM_2d
- **Issues**: https://github.com/HYGUO1993/FEM_2d/issues
- **文档**: 见仓库中的 README.md

---

**祝使用愉快！** 🚀
