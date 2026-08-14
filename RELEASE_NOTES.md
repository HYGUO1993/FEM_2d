# FemLab Studio Release v0.3.1

## 概述

**FemLab Studio** 是二维杆系/桁架/刚架有限元分析教学软件，包含：
- **Tauri 桌面应用**（`gui_tauri/`）— 交互式画布建模、求解、可视化一体
- **femcli**（`femcli.cpp`）— JSON 模型命令行求解器
- **FEM 核心引擎**（`barsystem.cpp`）— C++11 平面刚架/桁架求解

## 主要功能 (v0.3.1)

- ✅ 画布交互式建模：节点 / 杆件 / 荷载 / 约束 编辑
- ✅ **节点 CAD 式坐标输入**（节点工具下输入 `x,y` 或 `1000,500mm` 回车创建）
- ✅ **杆上集中力快捷交互**：点杆件上具体位置自动预选类型并填位置
- ✅ 求解：位移、杆端内力（N/V/M）、支座反力
- ✅ **精确内力图**：均布荷载下剪力线性变化、弯矩抛物线（含反弯点）；变形图 Hermite 曲线
- ✅ 显示选项卡：原图 / 轴力 N / 剪力 V / 弯矩 M / **挠度**
- ✅ LLM 自然语言建模助手（工具调用式 Agent，网格公式+模板提示，JSON 容错解析）
- ✅ 孤立节点诊断（求解失败原因明确提示）
- ✅ 10 个内置教程工程（含 **5x5 刚框架** 一键加载）
- ✅ 项目本地保存 / 导入导出 JSON 文件 / 撤销重做 / 中英双语 / 多主题

## 快速开始

### 桌面应用（推荐）

```powershell
# 开发运行
cd gui_tauri
npm install
npm run build
cd src-tauri
cargo tauri dev

# 或直接运行已构建的 exe
.\gui_tauri\src-tauri\target\debug\femlab-studio.exe
```

### 命令行求解 (femcli)

```bash
# 校验模型
build/bin/femcli.exe validate examples/beam_simple.json

# 求解
build/bin/femcli.exe solve examples/beam_simple.json -o result.json
```

### 构建 (Windows / Linux / macOS)

```bash
cmake -S . -B build
cmake --build build --parallel
ctest --test-dir build
```

## 模型格式

使用 JSON 模型（见 `docs/FEM_MODEL_SCHEMA.md` 与 `examples/`）：

```json
{
  "schemaVersion": "1.0",
  "nodes": [{ "id": 0, "type": "frame", "x": 0, "y": 0 }],
  "constraints": [{ "node": 0, "dofs": ["ux", "uy"] }],
  "elements": [{ "id": 0, "type": "frame", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 }],
  "materials": [{ "id": 0, "E": 210000000000, "mu": 0.3, "alpha": 0 }],
  "sections": [{ "id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.1 }],
  "loads": [{ "type": "nodalForce", "direction": "y", "value": -10000, "node": 1 }]
}
```

## 目录结构

```
FEM_2d/
├── gui_tauri/          # Tauri 桌面应用 (React 19 + Vite 8 + zustand 5)
├── barsystem.cpp/h     # FEM 核心求解引擎
├── femcli.cpp          # JSON 模型命令行求解器
├── CMakeLists.txt      # 构建配置
├── docs/               # 文档 (schema / LLM 集成 / 教程)
├── examples/           # JSON 模型示例
├── tests/              # 单元测试
├── verify/             # 求解器验证脚本
└── third_party/        # 第三方头文件 (nlohmann/json)
```

## 技术栈

| 层 | 技术 |
|---|---|
| 前端 | React 19, Vite 8, zustand 5, SVG 画布 |
| 桌面壳 | Tauri v2 (WebView2) |
| 后端 | Rust (IPC), C++11 (FEM 核心) |
| 构建 | CMake, cargo, npm |

## 测试

```bash
# 单元测试
ctest --test-dir build

# 求解器验证 (femcli vs 独立参考解, 13 案例)
python verify/compare_femcli_ref.py

# 渲染几何测试
node gui_tauri/scripts-test/render_check.mjs
```

## 版本历史

### v0.3.1 (2026-08-14)
- **求解器修复**：momentOnPoint 节点弯矩崩溃(堆损坏)、轴向力/均布被忽略、多荷载 FEF 覆盖、truss 节点 rz 垃圾值
- **内力图重做**：精确绘制杆内分布（均布荷载剪力线性/弯矩抛物线含反弯点），修复"缩成竖线"缺陷；变形图 Hermite 三次插值
- **显示选项卡**：新增「挠度」，变形图不再常驻画布
- **建模**：节点坐标输入(CAD式)、杆上集中力点击定位、5x5 刚框架教程、孤立节点诊断
- **LLM Agent**：修复多轮工具调用 400(Rust 消息字段丢失)、网格公式+模板提示、JSON 容错解析、建模成功率实测提升
- **构建**：femcli 同步副本自动化(CMake POST_BUILD)、golden 参考求解器 LDLT bug 修复

### v0.3.0 (2026-08-13)
- LLM Agent 工具调用架构 (建模/求解/解读闭环)
- /help 命令 + 9 个内置教程工程一键加载
- 全界面中英文双语切换 + 多主题 (深色/浅色/海洋)
- 修复: LDLT 回代 bug、均布荷载案例黑屏、校验误报、导出PNG原生对话框

### v0.2.0 (2026-08-13)
- 全新 Tauri 桌面端：画布建模 / 求解 / N·V·M 内力图 / LLM 辅助 / 撤销重做
- femcli JSON 命令行求解器
- 旧 PyQt6 GUI 退役归档（`_archive_local/`）

### v1.0.1 (2026-04-18)
- 修复 CI 构建失败，添加 .gitignore

### v1.0 (2026-04-18)
- C++ 核心引擎 + Python 可视化 + 基础 GUI

## 许可证

请参见项目仓库的 LICENSE 文件

## 联系方式

- GitHub: https://github.com/HYGUO1993/femlab-studio
- Issues: https://github.com/HYGUO1993/femlab-studio/issues
