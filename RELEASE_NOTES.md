# FemLab Studio Release v0.3.0

## 概述

**FemLab Studio** 是二维杆系/桁架/刚架有限元分析教学软件，包含：
- **Tauri 桌面应用**（`gui_tauri/`）— 交互式画布建模、求解、可视化一体
- **femcli**（`femcli.cpp`）— JSON 模型命令行求解器
- **FEM 核心引擎**（`barsystem.cpp`）— C++11 平面刚架/桁架求解

## 主要功能 (v0.3.0)

- ✅ 画布交互式建模：节点 / 杆件 / 荷载 / 约束 编辑
- ✅ 求解：位移、杆端内力（N/V/M）、支座反力
- ✅ 结果可视化：变形图 + **轴力图 / 剪力图 / 弯矩图** 叠加显示
- ✅ 节点位移数值标注（挠度直观显示）
- ✅ LLM 自然语言建模助手（OpenAI 兼容 API）
- ✅ 项目本地保存 / 导入导出 JSON 文件
- ✅ 撤销 / 重做（Ctrl+Z / Ctrl+Y）
- ✅ 快捷键：Delete 删除选中、Ctrl+Enter 求解
- ✅ 可自定义布局：三栏拖拽调宽、LLM 面板右侧/底部切换
- ✅ 完整安装包：MSI + NSIS setup.exe

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

# 求解器验证
python verify/verify_femcli.py
```

## 版本历史

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
