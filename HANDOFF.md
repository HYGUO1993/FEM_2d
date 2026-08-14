# FemLab Studio — 项目接手文档 (HANDOFF)

> 写给下一个开发者 / AI Agent。读完本文档 + 跑通构建，即可接手。

## 0. 项目速览

**FemLab Studio**：二维杆系有限元分析（平面刚架/桁架）教学软件。
- **前端**：React 19 + Vite 8 + zustand 5（`gui_tauri/src/`）
- **桌面壳**：Tauri v2（`gui_tauri/src-tauri/`）
- **求解核心**：C++11（`barsystem.cpp` / `femcli.cpp`）
- **LLM 助手**：OpenAI 兼容 API（deepseek 等），工具调用式 Agent

当前版本：**v0.3.0**（GitHub Release 已发布，NSIS 安装包）
仓库：https://github.com/HYGUO1993/femlab-studio

---

## 1. 架构总览

```
┌────────────────────────────────────────────────┐
│ React 前端 (gui_tauri/src/)                     │
│  App.jsx       布局/顶栏/设置面板/快捷键         │
│  store.js      zustand 全局状态(模型/历史/偏好)  │
│  agent.js      LLM Agent 工具调用循环            │
│  canvas/       SVG 渲染 (render/transform)      │
│  components/   Sidebar/CanvasView/Inspector/    │
│                ResultsPanel/ChatPanel/Toolbar   │
│  i18n.js       zh/en 双语字典                    │
│  model.js      FEM 模型 schema + 校验            │
│  ipc.js        Tauri IPC 封装                    │
└──────────┬─────────────────────────────────────┘
           │ invoke
┌──────────▼─────────────────────────────────────┐
│ Rust 后端 (gui_tauri/src-tauri/src/lib.rs)      │
│  solve_model     调 femcli.exe 求解             │
│  llm_chat_tools  LLM 对话(支持 tools)           │
│  save/load/list  项目管理                        │
│  export/import   文件导入导出(rfd 对话框)        │
└──────────┬─────────────────────────────────────┘
           │ spawn 子进程
┌──────────▼─────────────────────────────────────┐
│ C++ 求解器 (femcli.cpp + barsystem.cpp)         │
│  JSON 模型 → 有限元求解 → JSON 结果              │
│  (skyline 刚度 / LDLT / 固端力法)               │
└────────────────────────────────────────────────┘
```

---

## 2. 构建与运行（Windows）

### 2.1 依赖
- **Node.js ≥ 20**（开发用 25.8）
- **Rust**（stable，tauri-cli 2.x）
- **MSVC**（VS 2022，C++ 工具链）
- **CMake ≥ 3.10**
- Python（可选，辅助脚本）

### 2.2 完整构建流程

```powershell
# 1) C++ 求解器 (femcli.exe)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
# 产物: build/bin/Release/femcli.exe
#        + build/bin/femcli.exe (CMake POST_BUILD 自动同步, tauri 打包引用此副本)

# 2) 前端
cd gui_tauri
npm install
npm run build          # 产物: gui_tauri/dist/

# 3) Tauri 打包
cd src-tauri
$env:WIX = "$env:LOCALAPPDATA\tauri\WixTools311"   # 固定 WiX 3.11
cargo tauri build      # 正式安装包 (NSIS)
# 或
cargo tauri build --debug   # 调试 exe (快速迭代)
```

### 2.3 运行
- 调试版：`gui_tauri/src-tauri/target/debug/femlab-studio.exe`
- 正式安装：NSIS 安装包 → `%LOCALAPPDATA%\FemLab Studio\femlab-studio.exe`

### 2.4 关键路径
| 用途 | 路径 |
|---|---|
| femcli 求解器 | `build/bin/femcli.exe`（Tauri 打包为 resource） |
| 安装包产物 | `gui_tauri/src-tauri/target/release/bundle/nsis/` |
| 教学案例 | `examples/*.json`（9 个，前端内置副本 `gui_tauri/src/tutorialProjects.js`） |

---

## 3. 数据模型（JSON Schema）

```json
{
  "schemaVersion": "1.0",
  "title": "名称",
  "solver": "builtin",
  "nodes": [{ "id": 0, "type": "frame|truss", "x": 0.0, "y": 0.0 }],
  "constraints": [{ "node": 0, "dofs": ["ux","uy","rz"] }],
  "elements": [{ "id": 0, "type": "frame|truss", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 }],
  "materials": [{ "id": 0, "E": 210000000000, "mu": 0.3, "alpha": 0 }],
  "sections": [{ "id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.1 }],
  "loads": [{ "type": "...", "direction": "x|y|rz", "value": 0, "node": 0, "element": 0, "position": 0, "T0": 0, "T1": 0 }]
}
```

**荷载类型**（`femcli.cpp LoadTypeFromStr`）：
- `nodalForce` 节点集中力（node）
- `lateralForce` 杆件横向集中力（element+position）
- `lateralUniformPressure` 横向均布（element+position=作用长度）
- `lateralLinearlyPressure` 横向线性分布（element）
- `momentOnPoint` 节点弯矩（node；⚠️ 曾因求解器按 element 处理而崩溃，2026-08 已修）
- `axialForce` / `axialPressure` 轴向力/均布（element+position；⚠️ 曾被求解器静默忽略，2026-08 已实现 FEF）
- `temperature` 温度（T0/T1）
- `supportMove` 支座位移（node，仅受约束节点）

**结果 JSON**：`displacements`（ux/uy/rz）、`endForces`（nodeI/nodeJ 的 N/V/M 局部坐标）、`reactions`、`stats`。

---

## 4. 求解核心（C++）

- `barsystem.cpp/h`：有限元核心（组装/求解/后处理）
  - skyline 下三角刚度存储 + LDLT 求解
  - **已修复**：LDLT 回代曾用 `z[j]` 而非 `pB[j]`，导致多自由度结构数值错误（2026-08 修复）
  - 固端力法支持单元荷载
  - `SupportMoveAssembly`：支座位移等效荷载
- `femcli.cpp`：JSON CLI（`solve`/`validate`），**MSVC /MT 静态链接**（仅依赖 KERNEL32.dll）
- 测试：`ctest --test-dir build -C Release`（unit_tests）、`python verify/verify_femcli.py`（golden 回归，输入在 `tests/golden/inputs/`）
  - ⚠️ `verify/golden_generate.py` 是独立 Python 参考求解器，曾含与 C++ 相同的 LDLT 回代 bug（`z[j]`→`x[j]`，2026-08 已修并重新生成 golden.json）；改它后必须重跑 `python verify/golden_generate.py`
  - 求解器对照：`python verify/compare_femcli_ref.py <model.json>...`（femcli vs 独立参考，逐项 maxdiff）
  - 前端渲染测试：`node gui_tauri/scripts-test/render_check.mjs`（内力分布公式/反弯点/变形图几何）

### 5.1 内力图渲染（重要，勿回退）
- `gui_tauri/src/canvas/render.js` `computeForceDiagram` 按单元荷载**精确绘制杆内分布**：
  `V(x)=V_i+Σq·x`、`M(x)=M_i-V_i·x-Σq·x²/2`（j 端自动等于 -M_j，跨节点连续）
  ——均布荷载下剪力线性、弯矩抛物线（含反弯点），曾误画为"杆内恒定 V_i + 端值直线"。
- `computeDeformed` 用 Hermite 三次插值（节点 rz 转角）画真实弯曲变形，勿改回直线。
- 符号约定：弯矩 M_i 逆时针为正（CCW, y 向上），局部 x' 沿杆轴、y' 左侧法向；
  荷载 value 为局部带符号值（与 `FixedEndForceCalcu` 一致）。

---

## 5. LLM Agent（当前实现）

- 后端 `llm_chat_tools(config, messages, tools)`：转发 tools 到 OpenAI 兼容 API
- 前端 `agent.js`：client-side 工具循环（最多 8 轮）
  - 工具：`get_current_model` / `validate_model` / `solve` / `get_result_summary` / `apply_model`
  - `apply_model` 当前自动执行（无确认卡，未来可加）
- `llmSystemPrompt.js`：`buildAgentSystemPrompt`（Agent 版）+ 旧版 `buildSystemPrompt`
- 配置：ChatPanel 设置面板（base_url/api_key/model），存 `%APPDATA%\com.femlab.studio\llm_config.json`
- 设计文档：`docs/LLM_AGENT_ARCHITECTURE.md`（Phase 2/3 未实现）

---

## 6. 已知问题 / 注意事项

1. **MSI 安装包不可用**（1603）：Tauri 生成的 MSI 在此 Win11 环境 LaunchCondition 表读取失败。
   → **只产 NSIS**（`tauri.conf.json` `targets: ["nsis"]`）。不要改回 "all"。
2. **上传 GitHub Release 附件**：必须用**直接二进制 body**（`Content-Type: application/octet-stream`），
   multipart 上传会被污染（文件变成 multipart 报文）。用 `docs/scripts/upload_release.py` 思路。
3. **femcli 必须静态链接**（/MT）：CMakeLists.txt 已强制，勿改回 /MD，否则分发缺运行库。
4. **i18n 循环依赖陷阱**：`render.js`（纯几何）**禁止 import i18n**——用 labelKey/labelArgs 语义键，
   CanvasView 渲染时翻译。否则 ESM TDZ 崩溃（黑屏）。
5. **均布荷载案例**：`lateralUniformPressure` 需 `element` 有效且 `node=-1`；校验按类型分流。
6. **WebView2 调试**：window title 不反映 document.title；JS 错误用 Rust 写文件诊断
   （临时加 debug_write 命令，用完删除）。
7. **版本号三处同步**：`tauri.conf.json` / `Cargo.toml` / `package.json`。
8. **CI**：`ci.yml`（C++ 构建）绿；`release.yml` 触发 `fem-v*`（勿用 `v*`，与 Tauri 发布冲突）；
   `tauri-release.yml` 触发 `app-v*`（Tauri 自动发布，尚未实测）。
9. **桁架节点 rz**：truss 节点无转角自由度，`DOFIndexCalcu` 会把 `iaDOFIndex[2]` 置 -1
   （曾遗留 0 → femcli 结果 JSON 中 truss 节点 rz 输出自由度 0 的位移垃圾值，2026-08 已修）；
   改动 DOF 编号后必须重跑 `verify/compare_femcli_ref.py`（含 truss_bridge 等桁架案例）。

---

## 7. 最近改动（2026-08）

- v0.3.0：LLM Agent 工具循环、/help+9教程工程、全界面中英双语、三主题、显示选项
- 修复：LDLT 回代 bug（核心数值）、均布案例黑屏（ESM 循环）、校验误报、导出 PNG 原生对话框
- 打包：NSIS-only、femcli 静态链接、femcli 打包为 resource
- 性能：CanvasView useMemo 缓存、结果表格限 500 行

---

## 8. 下一步建议（按优先级）

1. **LLM Agent Phase 2**：结果解读（对比理论值 PL/qL²/8/5qL⁴/384EI）、apply_model 确认卡
2. **femcli 常驻服务**：stdin/stdout 管道通信，消除每次 spawn 开销（当前 ~10ms 可接受，非紧急）
3. **流式输出**：LLM 回复 SSE 打字机
4. **教学增强**：参数扫描动画、PDF 报告导出、型钢/材料库
5. **MSI 修复**：调查 Tauri MSI 在 Win11 的 LaunchCondition 问题（可暂缓，NSIS 够用）

---

## 9. 常用命令速查

```powershell
# 求解测试
build/bin/femcli.exe solve examples/simply_supported_udl.json -o result.json
build/bin/femcli.exe validate examples/beam_simple.json

# 前端快速迭代
cd gui_tauri && npm run dev          # Vite dev server (需 tauri dev 才接壳)
cd gui_tauri/src-tauri && cargo tauri dev

# 打包
cd gui_tauri/src-tauri
$env:WIX = "$env:LOCALAPPDATA\tauri\WixTools311"
cargo tauri build --debug    # 调试 exe
cargo tauri build            # NSIS 安装包

# 测试
ctest --test-dir build
python verify/verify_femcli.py
```

## 10. 联系方式 / 资源

- 仓库：https://github.com/HYGUO1993/femlab-studio
- Release：https://github.com/HYGUO1993/femlab-studio/releases
- 架构设计：`docs/LLM_AGENT_ARCHITECTURE.md`、`docs/ROADMAP.md`
- 模型格式：`docs/FEM_MODEL_SCHEMA.md`
