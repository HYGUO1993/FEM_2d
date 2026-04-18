# FEM_2d 开发文档

## 阶段一进展报告：基础架构

### 已完成的工作

#### 1. C++ 代码模块化 ✅

已将原始的单文件结构重构为模块化架构：

```
src/
├── fem_types.h      # 数据结构定义（Material, Section, Node, Element等）
├── fem_core.h       # Model类声明（面向对象封装）
├── fem_core.cpp     # Model类实现
├── fem_api.h        # C API接口声明
└── fem_api.cpp      # C API实现
```

**特点：**
- 使用 `FEM` 命名空间避免命名冲突
- `Model` 类封装所有FEM逻辑
- RAII 内存管理，自动清理资源
- 完整的错误处理机制

#### 2. C API 接口层 ✅

提供完整的 C 接口用于跨语言互操作：

```c
// 模型管理
FEMModelHandle FEM_CreateModel();
void FEM_DestroyModel(FEMModelHandle model);

// 模型构建
int FEM_AddNode(FEMModelHandle model, int type, double x, double y);
int FEM_AddElement(...);
int FEM_AddMaterial(...);
// ... 等

// 分析
int FEM_Solve(FEMModelHandle model);

// 结果获取
const double* FEM_GetDisplacements(FEMModelHandle model, int* count);
const double* FEM_GetForces(...);
const double* FEM_GetReactions(...);
```

**特点：**
- 不透明句柄（opaque handle）设计
- 完整的错误信息反馈
- 内存安全，防止泄漏

#### 3. Python 绑定 ✅

创建了两种 Python 接口：

**方案 A：pybind11 绑定（推荐）**
```python
# python/pyfem_bindings.cpp
import pyfem

model = pyfem.Model()
model.add_material(E=210e9, mu=0.3, alpha=1.2e-5)
model.add_node(pyfem.FRAME_NODE, 0.0, 0.0)
model.solve()
displacements = model.get_displacements()  # 返回 NumPy array
```

**方案 B：ctypes 封装（备选）**
```python
# python/fem2d/model_ctypes.py
from fem2d import Model, FRAME_NODE

model = Model()  # 自动选择可用的后端
model.add_material(E=210e9, mu=0.3, alpha=1.2e-5)
```

**特点：**
- 自动类型转换
- NumPy 数组集成
- Pythonic API 设计
- 完整的文档字符串

#### 4. 构建系统升级 ✅

更新了 `CMakeLists.txt`：

```cmake
# 选项
option(BUILD_SHARED_LIB "Build FEM as a shared library" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

# 生成共享库
add_library(fem_core SHARED ${FEM_CORE_SOURCES})

# Python 绑定（可选）
if(BUILD_PYTHON_BINDINGS)
    pybind11_add_module(pyfem python/pyfem_bindings.cpp)
endif()
```

**特点：**
- 支持构建共享库（.so/.dll）
- 可选的 Python 绑定编译
- 保持向后兼容（原始 fem_run 仍可用）

#### 5. Python 包结构 ✅

```
python/
├── fem2d/
│   ├── __init__.py          # 包入口，自动选择后端
│   └── model_ctypes.py      # ctypes 封装
├── gui/
│   ├── __init__.py
│   └── main_window.py       # PySide6 主窗口
├── fem_gui.py               # GUI 应用入口
├── test_api.py              # API 测试脚本
└── pyfem_bindings.cpp       # pybind11 绑定源码
```

#### 6. GUI 框架基础 ✅

创建了 PySide6 GUI 应用框架：

```python
# 运行 GUI
python python/fem_gui.py
```

**已实现的功能：**
- 主窗口布局（菜单栏、工具栏、状态栏）
- 模型树面板（左侧，显示节点/单元/载荷等）
- 属性面板（右侧，显示选中对象属性）
- 中央可视化区域（占位符，待集成 PyVista）
- 文件操作菜单（新建、打开、保存）
- 分析菜单（求解按钮）

### 构建与测试

#### 构建 C++ 库

```bash
# Linux/macOS
mkdir build && cd build
cmake .. -DBUILD_SHARED_LIB=ON
make

# 可选：构建 Python 绑定
cmake .. -DBUILD_SHARED_LIB=ON -DBUILD_PYTHON_BINDINGS=ON
make
```

#### 测试 Python API

```bash
# 确保库已构建
cd python
python test_api.py

# 预期输出：
# ✓ Successfully imported fem2d
# ✓ Model created
# ✓ Material added
# ... (更多测试)
```

#### 运行 GUI（预览）

```bash
# 安装依赖
pip install PySide6 pyvista matplotlib numpy

# 运行
python python/fem_gui.py
```

### 技术架构图

```
┌─────────────────────────────────────────────────────────┐
│                    Python GUI (PySide6)                 │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  │
│  │  Main Window │  │  Visualizer  │  │   Dialogs    │  │
│  └──────────────┘  └──────────────┘  └──────────────┘  │
└───────────────────────┬─────────────────────────────────┘
                        │
┌───────────────────────┴─────────────────────────────────┐
│              Python Bindings Layer                      │
│   ┌─────────────────┐        ┌──────────────────┐      │
│   │  pybind11 API   │   OR   │   ctypes API     │      │
│   └─────────────────┘        └──────────────────┘      │
└───────────────────────┬─────────────────────────────────┘
                        │
┌───────────────────────┴─────────────────────────────────┐
│                   C API Layer (fem_api.h)               │
│   FEM_CreateModel(), FEM_Solve(), FEM_GetResults()...  │
└───────────────────────┬─────────────────────────────────┘
                        │
┌───────────────────────┴─────────────────────────────────┐
│            C++ Core Library (fem_core.cpp)              │
│                   class FEM::Model                      │
│   ┌──────────┐ ┌──────────┐ ┌──────────┐              │
│   │ Geometry │ │  Solver  │ │   I/O    │              │
│   └──────────┘ └──────────┘ └──────────┘              │
└─────────────────────────────────────────────────────────┘
                        │
┌───────────────────────┴─────────────────────────────────┐
│          Legacy Code (barsystem.cpp) - Reused           │
│   函数：LDLTSolve, GKAssembly, InternalForceCalcu...   │
└─────────────────────────────────────────────────────────┘
```

### 下一步计划

#### 阶段一剩余工作：
- [ ] 修复编译错误（src/fem_core.cpp 中的包含路径）
- [ ] 实现 `Model::saveToFile()` 方法
- [ ] 添加单元测试覆盖新 API
- [ ] 编写 API 使用文档

#### 阶段二：GUI 增强
- [ ] 集成 PyVista 3D 可视化
- [ ] 实现模型树的交互（点击显示属性）
- [ ] 添加工具面板（添加节点/单元的工具）
- [ ] 实现文件 I/O（读取 .txt 格式）

### 文件清单

**新增文件：**
- `src/fem_types.h` (150 行) - 数据结构
- `src/fem_core.h` (105 行) - Model 类声明
- `src/fem_core.cpp` (300 行) - Model 类实现
- `src/fem_api.h` (200 行) - C API 声明
- `src/fem_api.cpp` (300 行) - C API 实现
- `python/pyfem_bindings.cpp` (250 行) - pybind11 绑定
- `python/fem2d/__init__.py` (60 行) - Python 包
- `python/fem2d/model_ctypes.py` (200 行) - ctypes 封装
- `python/gui/main_window.py` (350 行) - GUI 主窗口
- `python/gui/__init__.py` (10 行)
- `python/fem_gui.py` (80 行) - GUI 入口
- `python/test_api.py` (200 行) - 测试脚本
- `DEVELOPMENT.md` (本文件)

**修改文件：**
- `CMakeLists.txt` - 添加共享库和 Python 绑定支持

**总计新增代码：** 约 2200 行

### 验收标准

#### 已完成 ✅：
- [x] C++ 代码模块化（命名空间、类封装）
- [x] C API 接口完整定义
- [x] Python 绑定框架（pybind11）
- [x] Python ctypes 封装
- [x] CMake 构建系统更新
- [x] GUI 框架基础（PySide6）
- [x] 测试脚本

#### 待完成 ⏳：
- [ ] 编译通过并生成共享库
- [ ] Python API 测试通过
- [ ] GUI 可以启动并显示窗口

### 技术决策

#### 为何选择 pybind11？
1. **类型安全**：自动类型转换，编译期检查
2. **性能**：零开销抽象，接近原生 C++ 速度
3. **易用性**：简洁的绑定语法，自动文档生成
4. **NumPy 集成**：原生支持 NumPy 数组
5. **现代 C++**：利用 C++11/14 特性

#### 为何同时提供 ctypes 封装？
1. **备选方案**：不依赖 pybind11 编译
2. **简单部署**：只需共享库，无需 Python 扩展
3. **调试方便**：纯 Python 代码易于修改

#### 为何选择 PySide6？
1. **官方支持**：Qt 官方 Python 绑定
2. **LGPL 协议**：适合开源项目
3. **功能完整**：完整的 Qt 功能支持
4. **跨平台**：Windows/Linux/macOS 一致体验

### 已知问题

1. **编译问题**：`src/fem_core.cpp` 需要调整包含路径
   - 应该包含 `"../barsystem.h"` 而不是 `"barsystem.h"`

2. **链接问题**：共享库需要正确链接所有函数
   - 可能需要将 `barsystem.cpp` 的函数拆分或导出

3. **Python 路径**：ctypes 查找库的路径可能需要调整
   - 可通过环境变量 `LD_LIBRARY_PATH` 或安装到系统路径解决

### 开发环境

**推荐配置：**
- **OS**: Ubuntu 20.04+ / macOS 12+ / Windows 10+
- **编译器**: GCC 9+ / Clang 10+ / MSVC 2019+
- **CMake**: 3.10+
- **Python**: 3.8+
- **依赖库**:
  - PySide6 6.0+
  - PyVista 0.38+
  - Matplotlib 3.5+
  - NumPy 1.21+
  - pybind11 2.10+ (可选)

### 贡献指南

欢迎贡献！请遵循以下步骤：

1. Fork 仓库
2. 创建特性分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 开启 Pull Request

**编码规范：**
- C++: Google C++ Style Guide
- Python: PEP 8
- 注释: Doxygen (C++) / NumPy docstring (Python)

---

**文档版本**: 1.0
**最后更新**: 2026-04-18
**作者**: FEM_2d 开发团队
