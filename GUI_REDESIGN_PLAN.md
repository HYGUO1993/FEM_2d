# FemLab Studio GUI 重构实施计划

> 目标执行者：按本计划零来回实现。所有决策已定死，不要重新设计，照做即可。
> 技术栈已由用户拍板：**React 18 + Vite 5 + zustand**（前端）；**Tauri v2 + ureq**（后端）。
> 范围：**P1（画布交互）+ P2（项目持久化）+ P3（LLM 自然语言建模）全做**。

---

## 0. 目标

把现有「JSON 编辑器 + 结果表格」两栏界面，重构为类 opencode/codex 的工程化界面：

- **左栏**：项目列表（多项目新建/切换/删除，本地 JSON 持久化）
- **中央**：SVG 画布（手动画节点/杆件/荷载/约束，pan/zoom，求解后叠加变形图）
- **右栏**：工具栏（交互模式切换）+ 属性面板（选中项）+ 结果面板 + LLM 对话面板
- **LLM**：OpenAI 兼容 API，自然语言建模 → 解析成模型 → 上画布 → 自动求解

---

## 1. 现状（务必先理解，别破坏现有可用的求解链路）

- 前端现在在 `gui_tauri/src/`（纯静态 `index.html` + `main.js` + `styles.css`），会被本次重构**整体替换**。
- 后端 `gui_tauri/src-tauri/src/lib.rs` 现有两个命令，**保留不动**：
  - `solve_model(model_json: String) -> Result<String, String>`：写临时文件 → 调 `femcli.exe solve` → 返回结果 JSON 字符串。
  - `ping() -> String`。
- 求解内核 `femcli.exe` 在 `build/bin/`，`locate_femcli` 会在 exe 同目录（`target/debug/`）找到它，不要动这个逻辑。
- `tauri.conf.json` 当前 `frontendDist: "../src"`，`withGlobalTauri: true`，`custom-protocol` 是默认 feature（保证 `cargo build` 出的 exe 内嵌前端资源，无需 dev server）。

## 2. FEM 模型 schema（前后端唯一数据契约，画布与 femcli 都遵循它）

内存模型对象 = 可直接 `JSON.stringify` 后传给 `solve_model` 的结构，**字段名必须完全一致**：

```json
{
  "schemaVersion": "1.0",
  "title": "字符串标题",
  "solver": "builtin",
  "nodes":       [ { "id": 0, "type": "frame|truss", "x": 0.0, "y": 0.0 } ],
  "constraints": [ { "node": 0, "dofs": ["ux","uy","rz"] } ],
  "elements":    [ { "id": 0, "type": "frame|truss", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 } ],
  "materials":   [ { "id": 0, "E": 210000000000.0, "mu": 0.3, "alpha": 0.0 } ],
  "sections":    [ { "id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.1 } ],
  "loads":       [ { "type": "nodalForce", "direction": "x|y|rz", "value": -10000.0, "node": 1 } ]
}
```

约束语义：`dofs` 里出现的分量被约束（位移=0）。荷载方向 `x/y/rz` 对应节点自由度。

femcli 结果 JSON 字段（前端结果面板要读这些，字段名已确认）：

```json
{
  "status": "ok|error",
  "message": "…",
  "stats": { "nodeCount": 3, "elementCount": 2, "freeDOF": 6, "totalDOF": 9 },
  "displacements": [ { "node": 0, "ux": 0.0, "uy": 0.0, "rz": 0.0 } ],
  "endForces":    [ { "element": 0, "nodeI": {"N":0,"V":0,"M":0}, "nodeJ": {"N":0,"V":0,"M":0} } ],
  "reactions":    [ { "node": 0, "ux": 0.0, "uy": 0.0, "rz": 0.0 } ]
}
```

---

## 3. 后端改动（`gui_tauri/src-tauri/`）

### 3.1 `Cargo.toml` 加依赖

在 `[dependencies]` 里追加（保持现有 tauri/serde/serde_json 不变）：

```toml
ureq = { version = "2", features = ["json"] }
```

（用 ureq 而非 reqwest：纯阻塞、体积小、不引入 tokio 重依赖。）

### 3.2 `lib.rs` 新增命令（追加，不动 `solve_model`/`ping`/`locate_femcli`）

新增以下 `#[tauri::command]`，并全部加入 `generate_handler!`：

```rust
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Default)]
struct LlmConfig {
    #[serde(default)] base_url: String,   // 例 "https://api.deepseek.com/v1"
    #[serde(default)] api_key: String,
    #[serde(default)] model: String,      // 例 "deepseek-chat"
}

#[derive(Serialize, Deserialize)]
struct ChatMessage {
    role: String,     // "system" | "user" | "assistant"
    content: String,
}
```

命令签名与实现要点：

| 命令 | 签名 | 行为 |
|---|---|---|
| `llm_chat` | `(config: LlmConfig, messages: Vec<ChatMessage>) -> Result<String, String>` | POST `{base_url}/chat/completions`，body `{"model","messages":[...],"temperature":0.2,"stream":false}`，header `Authorization: Bearer {api_key}`。返回 `choices[0].message.content`。网络/解析错误一律 `Err(字符串说明)`。 |
| `save_project` | `(name: String, model_json: String) -> Result<(), String>` | 写 `<app_data_dir>/projects/<name>.json`（先建目录）。 |
| `list_projects` | `() -> Result<Vec<String>, String>` | 列出 `<app_data_dir>/projects/` 下 `*.json` 文件名（去 `.json` 后缀，按名排序）。目录不存在返回空数组。 |
| `load_project` | `(name: String) -> Result<String, String>` | 读 `<app_data_dir>/projects/<name>.json` 返回原文。 |
| `delete_project` | `(name: String) -> Result<(), String>` | 删对应文件；不存在则 Ok（幂等）。 |
| `get_llm_config` | `() -> Result<LlmConfig, String>` | 读 `<app_data_dir>/llm_config.json`；不存在返回默认（空）LlmConfig。 |
| `set_llm_config` | `(config: LlmConfig) -> Result<(), String>` | 写 `<app_data_dir>/llm_config.json`（先建目录）。 |

`app_data_dir` 获取：`app.path().app_data_dir()`（需要 `app: tauri::AppHandle` 参数，仿照 `locate_femcli` 的写法；没有 app 参数的命令用 `tauri::Manager` trait）。文件名里 `<name>` 需做基础清洗（替换 `\ / : * ? " < > |` 为 `_`），防止路径穿越。

`generate_handler!` 最终应为：

```rust
.invoke_handler(tauri::generate_handler![
    solve_model, ping,
    llm_chat, save_project, list_projects, load_project, delete_project,
    get_llm_config, set_llm_config
])
```

### 3.3 验证（后端独立验证）

`cargo check` 必须通过。**注意**：`cargo check` 会跑 build script，`frontendDist` 改为 `../dist` 后需要 `dist` 目录存在。若 `gui_tauri/dist/` 尚不存在，先手动 `mkdir dist` 并放一个空 `dist/index.html` 占位，再 `cargo check`。

---

## 4. 前端改动（`gui_tauri/`，React + Vite）

### 4.1 工程文件

**`gui_tauri/package.json`**：

```json
{
  "name": "femlab-studio-frontend",
  "private": true,
  "version": "0.1.0",
  "type": "module",
  "scripts": {
    "dev": "vite",
    "build": "vite build",
    "preview": "vite preview"
  },
  "dependencies": {
    "react": "^18.3.1",
    "react-dom": "^18.3.1",
    "zustand": "^4.5.5"
  },
  "devDependencies": {
    "@vitejs/plugin-react": "^4.3.4",
    "vite": "^5.4.11"
  }
}
```

**`gui_tauri/vite.config.js`**：

```js
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  server: { port: 5173, strictPort: true },
  build: { outDir: "dist", target: "chrome105", sourcemap: false },
});
```

**`gui_tauri/index.html`**（Vite 入口，放根目录）：

```html
<!DOCTYPE html>
<html lang="zh-CN">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>FemLab Studio</title>
  </head>
  <body>
    <div id="root"></div>
    <script type="module" src="/src/main.jsx"></script>
  </body>
</html>
```

**`gui_tauri/.gitignore`**：追加 `node_modules`、`dist`。

**删除旧文件**：`gui_tauri/src/index.html`、`gui_tauri/src/main.js`、`gui_tauri/src/styles.css`（旧的纯静态前端整体作废，`src/` 改为放 React 源码）。

### 4.2 源码结构

```
gui_tauri/src/
├── main.jsx            # ReactDOM.createRoot 挂载 <App/>
├── App.jsx             # 三栏布局骨架
├── store.js            # zustand store（唯一状态源）
├── model.js            # 默认示例模型 + 工厂函数
├── ipc.js              # tauri invoke 封装
├── components/
│   ├── Sidebar.jsx     # 左栏：项目列表
│   ├── CanvasView.jsx  # 中央：SVG 画布
│   ├── Toolbar.jsx     # 右栏顶部：交互模式按钮
│   ├── Inspector.jsx   # 右栏：选中项属性编辑
│   ├── ResultsPanel.jsx# 右栏：求解结果表格
│   └── ChatPanel.jsx   # 右栏：LLM 对话 + 配置
├── canvas/
│   ├── transform.js    # 坐标变换 + pan/zoom
│   └── render.js       # 模型 → SVG 元素几何计算（含变形图）
└── styles.css          # 深色主题
```

### 4.3 `store.js`（zustand，状态与动作定义）

状态字段：`projects: string[]`、`currentProject: string`、`model`（schema 对象，见 §2）、`tool: 'select'|'node'|'element'|'load'|'constraint'|'erase'`、`selection: {type:'node'|'element'|'load'|'constraint', id:number} | null`、`results: object|null`、`solved: boolean`、`pendingNodeA: number|null`（element 模式第一个点）、`llmConfig: {base_url,api_key,model}`、`chatMessages: {role,content}[]`、`deformScale: number`。

动作（actions）：`newProject(name)`、`switchProject(name)`、`deleteProject(name)`、`renameProject`（可选）、`setTool`、`select(sel)`、`addNode(x,y)`、`moveNode(id,x,y)`、`deleteNode(id)`、`addElement(nodeI,nodeJ)`、`deleteElement(id)`、`addLoad(node,direction,value)`、`deleteLoad(id)`、`toggleConstraint(node,dof)`、`updateMaterial`/`updateSection`（属性面板编辑）、`setResults`、`setLlmConfig`、`pushChat`、`setDeformScale`。

关键规则（在 action 里强制）：
- 删除节点时，级联删除引用它的 elements / loads / constraints。
- 新增模型至少含 1 个默认 material（E=2.1e11, mu=0.3, alpha=0）和 1 个默认 section（A=0.01, Iz=1e-5, height=0.1），id 均从 0 起；新节点默认 `type:"frame"`，id 自增取当前最大 +1。
- 求解前做前端校验：至少 2 节点、1 单元、有约束、无悬空引用；不通过则 toast 提示并阻止求解。

### 4.4 `model.js`

导出 `defaultModel()`：一个 2m 简支梁（3 节点 frame：x=0/1/2，约束 node0=[ux,uy]、node2=[uy]，2 根 frame 单元，跨中 node1 竖向 -10000 集中力），与现有 `main.js` 里的 EXAMPLE 一致，作为画布初始内容。

### 4.5 `ipc.js`

封装 `invoke(cmd, args)`，用 `window.__TAURI__.core.invoke`。另封装高层方法：`solve(model)`、`llmChat(config, messages)`、`saveProject(name, model)`、`listProjects()`、`loadProject(name)`、`deleteProject(name)`、`getLlmConfig()`、`setLlmConfig(config)`。

### 4.6 画布（`CanvasView.jsx` + `canvas/`）—— 核心

**用 SVG 而非 Canvas**（React 天然契合、事件好绑、矢量清晰）。

- 背景：深色网格（CSS `background-image: linear-gradient` 或 SVG `<pattern>`，网格间距随 zoom 自适应，可选简化为固定 20px）。
- 坐标变换：模型坐标 y 向上，屏幕 y 向下 → `worldToScreen(x,y)` 需翻转 y。状态 `{scale, panX, panY}`（存 store 或 CanvasView 本地 state）。
  - `sx = x*scale + panX`，`sy = -y*scale + panY`（panX/panY 为屏幕原点偏移）。
  - 初始 scale 自适应：根据模型包围盒，让结构居中并占画布约 60%（首屏及「重置视图」按钮触发）。
  - 滚轮缩放（以鼠标位置为锚点）、中键或空格拖拽平移、「重置视图」按钮。
- 渲染元素：
  - 节点：`<circle>` + 编号 `<text>`。truss 节点与 frame 节点外观可区分（frame 实心、truss 空心）。选中高亮（描边 accent 色）。
  - 杆件：`<line>`。truss 用虚线、frame 用实线（2px）。选中高亮。
  - 荷载：从节点出发的箭头（`<line>` + 箭头 `<polygon>`），方向 x/y/rz（rz 画成绕节点的小圆弧箭头），旁标数值文本。颜色区分（如橙色）。
  - 约束：固定端画实心三角贴节点，滚轴（单个 dof）画三角+圆，双 dof 约束画带斜线阴影。用简单符号，颜色（如蓝色）。
  - 变形图：`results` 存在且 `solved` 时，用 `displacements` 计算 `newPos = pos + u*deformScale`，把变形后的杆件用**红色虚线**叠加在原结构上；`deformScale` 由结果位移量级自动估算（`scale = 0.15*span/max|u|`，`max|u|` 为 0 则 0），右栏有滑块微调。
- 交互（依据 `tool`）：
  - `select`：点击命中节点/杆件 → 选中并在 Inspector 显示属性；拖拽节点改其 x/y（实时更新 model 并重绘）；空白拖拽平移视图。
  - `node`：点击空白处，把屏幕坐标反算成模型坐标，新增节点（吸附到 0.1m 网格可选）。
  - `element`：点第一个节点（高亮记录 `pendingNodeA`）→ 点第二个节点 → 建杆件并清空 pending。点空白取消。
  - `load`：点节点 → 弹简单内联对话框（方向 x/y/rz + 数值）→ 确认后加荷载。或点节点后在 Inspector 输入。
  - `constraint`：点节点 → 在 Inspector 显示 ux/uy/rz 三个 checkbox，勾选即 toggle。
  - `erase`：点节点/杆件/荷载/约束 → 删除（节点级联）。
  - 命中检测：节点按屏幕距离 < 12px；杆件按点到线段距离 < 8px（换算到屏幕坐标）。
- 状态栏（画布底部或 App 底部）：`N 节点 · M 单元 · K DOF · 求解耗时 · ✔/✘ 模型状态`。

### 4.7 右栏面板

- **Toolbar**：6 个图标按钮（select/node/element/load/constraint/erase），active 高亮。用 SVG 小图标或文字缩写均可。
- **Inspector**：根据 `selection` 显示对应编辑控件：
  - 节点：x、y（number input，实时更新）、type（frame/truss 下拉）。
  - 杆件：nodeI、nodeJ、type、section、material（下拉）。
  - 荷载：direction、value、node。
  - 约束：node + ux/uy/rz checkbox。
  - 无选中：显示「模型」全局设置（materials/sections 表格编辑）。
- **ResultsPanel**：复用现有结果渲染逻辑（stats 卡片 + 位移/端力/反力三张表），但注意**表格用 React 直接 map 渲染**（别再犯 `body.push` 那种错误）。`status==="error"` 时显示错误框。结果要能驱动变形图（setResults）。
- **ChatPanel**：聊天流（用户/助手气泡）+ 输入框 + 发送按钮 + 「设置」按钮（弹出 base_url/api_key/model 表单，保存到 `setLlmConfig` 后端）。助手回复中的 JSON 自动提取并「应用到画布」按钮，点击后 `setModel` + 自动 `solve`。

### 4.8 LLM 流程

1. 前端组装 system prompt（固定文案，写死在前端 `llmSystemPrompt.js` 或 ChatPanel 内）：告诉 LLM「你是二维杆系有限元建模助手，只输出符合给定 schema 的严格 JSON 模型，不要输出解释，不要用 markdown 代码块包裹」。prompt 内附上 §2 的 schema 说明和字段语义。
2. `chatMessages` 附带当前模型 JSON（作为上下文，让 LLM 在现有模型上改，或全新生成）。
3. 调 `llmChat(config, messages)` → 后端 ureq 请求 → 返回文本。
4. 前端从返回文本提取 JSON（正则取第一个 `{` 到最后一个 `}`），`JSON.parse` 成功则存入 `pendingLlmModel`，UI 显示「已生成模型：N 节点 M 单元 [应用到画布]」；解析失败则当作普通文本显示在聊天流。

### 4.9 左栏 Sidebar（项目列表）

- 顶部「+ 新建项目」按钮 → 输入名字（默认 `项目 N`）→ 新建并切换。
- 项目列表：每项高亮当前项；点击切换（切换前自动保存当前项目）；右键或悬停显示删除按钮。
- 自动保存：任何 model 变更后 debounce 500ms 调 `saveProject(currentProject, JSON.stringify(model))`。启动时 `listProjects()` 拉列表；若有项目则加载最近/第一个，否则用 `defaultModel()` 建「未命名项目」。

### 4.10 样式 `styles.css`

深色工程主题（参考 codex/opencode）：背景 `#0f1115` 层级、面板 `#16181d`、边框 `#262a33`、文字 `#e6e8ec`、次要 `#9aa0ab`、accent 蓝 `#4da3ff`（或 `#3b82f6`）、错误 `#f87171`、成功 `#34d399`。等宽字体用系统 monospace。布局用 CSS Grid：`grid-template-columns: 240px 1fr 320px`，header 和状态栏 `flex-shrink:0`，中间画布 `min-width:0; min-height:0` 且 `overflow:hidden`。整体无滚动条溢出。

---

## 5. `tauri.conf.json` 改动

```json
"build": {
  "frontendDist": "../dist",
  "beforeDevCommand": "npm run dev",
  "devUrl": "http://localhost:5173"
}
```

其余（productName/identifier/withGlobalTauri/窗口）不变。`custom-protocol` feature 保持默认（build 时嵌入 dist）。

---

## 6. 构建与验证流程（最终集成必须跑通）

```bash
cd gui_tauri
npm install            # npmmirror 已配置
npm run build          # 产出 dist/

cd src-tauri
cargo build            # 嵌入 dist，产出 target/debug/femlab-studio.exe
```

验证清单：
1. `npm run build` 无报错，`gui_tauri/dist/index.html` 存在。
2. `cargo build` exit 0。
3. 运行 `femlab-studio.exe`：界面为深色三栏布局；画布初始显示简支梁示例。
4. 点「求解」：结果面板出位移/端力/反力，画布叠加红色变形图，无 `body.push` 类报错。
5. 新建项目、切换、删除、重启后项目仍在（持久化）。
6. 配置 LLM 后自然语言建模能生成模型并求解（无 key 时后端返回明确错误提示，不崩溃）。

---

## 7. 关键决策记录（不要推翻）

- 前端：React 18 + Vite 5 + zustand，SVG 画布（非 Canvas）。深色主题。
- 后端：新增 7 个命令（llm_chat/save_project/list_projects/load_project/delete_project/get_llm_config/set_llm_config），HTTP 用 ureq 阻塞式，项目与配置存 `app_data_dir`。
- 模型 schema 与 femcli 结果字段见 §2，是唯一契约，字段名必须精确匹配。
- `solve_model`/`ping`/`locate_femcli` 现有逻辑不动。
- 前端不发 HTTP（LLM 走后端代理避免 CORS），capabilities 保持 `core:default` 不变。
- 结果表格用 React 直接 map，不用手写 DOM 拼接。
