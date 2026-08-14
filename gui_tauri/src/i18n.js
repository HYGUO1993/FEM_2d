// 轻量 i18n: zh / en 字典 + t() 函数
// 用法: import { t, setLang } from "./i18n.js"
// 支持占位符: t("key", {name: "xx"}) → 替换 {name}

const dict = {
  // 顶栏
  "llm.toBottom": { zh: "LLM 移到底部", en: "LLM to Bottom" },
  "llm.toRight": { zh: "LLM 移到右侧", en: "LLM to Right" },
  "llmPosTip": { zh: "切换 LLM 对话位置", en: "Toggle LLM panel position" },
  "undo": { zh: "撤销", en: "Undo" },
  "redo": { zh: "重做", en: "Redo" },
  "resetView": { zh: "重置视图", en: "Reset View" },
  "solve": { zh: "求解", en: "Solve" },
  "solving": { zh: "求解中…", en: "Solving…" },
  "undoTip": { zh: "撤销 (Ctrl+Z)", en: "Undo (Ctrl+Z)" },
  "redoTip": { zh: "重做 (Ctrl+Y)", en: "Redo (Ctrl+Y)" },
  "app.subtitle": { zh: "二维杆系有限元 · 画布建模 · LLM 辅助", en: "2D Frame FEM · Canvas Modeling · LLM Assisted" },
  "status.nodes": { zh: "{n} 节点", en: "{n} nodes" },
  "status.elements": { zh: "{n} 单元", en: "{n} elements" },
  "status.dofs": { zh: "{n} DOF", en: "{n} DOF" },
  "status.solveTime": { zh: "求解耗时 {t}s", en: "Solve time {t}s" },
  "status.modelOk": { zh: "✔ 模型有效", en: "✔ Model OK" },
  "unnamed": { zh: "未命名项目", en: "Untitled" },
  "initFailed": { zh: "初始化失败: {msg}", en: "Init failed: {msg}" },

  // 侧栏
  "newProject": { zh: "+ 新建项目", en: "+ New Project" },
  "saveLocal": { zh: "保存到本地文件", en: "Save to File" },
  "importLocal": { zh: "导入本地文件", en: "Import File" },
  "noProjects": { zh: "暂无项目", en: "No projects" },
  "create": { zh: "创建", en: "Create" },
  "projectName": { zh: "项目名称", en: "Project name" },
  "rename": { zh: "改名", en: "Rename" },
  "duplicate": { zh: "复制", en: "Duplicate" },
  "delete": { zh: "删除项目", en: "Delete project" },
  "saveLocalTip": { zh: "将当前模型保存为本地 JSON 文件", en: "Save current model to a local JSON file" },
  "importLocalTip": { zh: "从本地 JSON 文件导入模型", en: "Import a model from a local JSON file" },
  "projectDefault": { zh: "项目 {n}", en: "Project {n}" },
  "toast.exported": { zh: "已导出: {path}", en: "Exported: {path}" },
  "toast.exportFailed": { zh: "导出失败: {msg}", en: "Export failed: {msg}" },
  "toast.canceledSave": { zh: "用户取消了保存", en: "Save canceled" },
  "toast.canceledPick": { zh: "用户取消了选择", en: "Selection canceled" },
  "toast.invalidJson": { zh: "文件内容不是有效的模型 JSON", en: "File content is not a valid model JSON" },
  "toast.missingNodes": { zh: "模型缺少 nodes 数组", en: "Model is missing nodes array" },
  "importedModel": { zh: "导入模型", en: "Imported Model" },
  "toast.imported": { zh: "已导入: {name}", en: "Imported: {name}" },
  "toast.importFailed": { zh: "导入失败: {msg}", en: "Import failed: {msg}" },
  "toast.created": { zh: "已创建项目: {name}", en: "Created project: {name}" },
  "toast.switchFailed": { zh: "切换项目失败: {msg}", en: "Switch project failed: {msg}" },
  "toast.deleted": { zh: "已删除项目: {name}", en: "Deleted project: {name}" },
  "toast.deleteFailed": { zh: "删除失败: {msg}", en: "Delete failed: {msg}" },
  "toast.duplicated": { zh: "已复制: {name}", en: "Duplicated: {name}" },
  "toast.duplicateFailed": { zh: "复制失败: {msg}", en: "Duplicate failed: {msg}" },
  "toast.renameExists": { zh: "同名项目已存在", en: "A project with the same name exists" },
  "toast.renamed": { zh: "已改名: {name}", en: "Renamed: {name}" },
  "toast.renameFailed": { zh: "改名失败: {msg}", en: "Rename failed: {msg}" },
  "dupeSuffix": { zh: " (副本)", en: " (copy)" },
  "dupeSuffixN": { zh: " (副本{n})", en: " (copy{n})" },

  // 工具栏
  "tool.select": { zh: "选择", en: "Select" },
  "tool.node": { zh: "节点", en: "Node" },
  "tool.element": { zh: "杆件", en: "Member" },
  "tool.load": { zh: "荷载", en: "Load" },
  "tool.constraint": { zh: "约束", en: "Support" },
  "tool.erase": { zh: "删除", en: "Erase" },

  // 画布提示
  "hint.select": { zh: "选择: 点击节点/杆件选中, 拖拽节点移动, 空白拖拽平移", en: "Select: click node/member, drag node to move, drag empty area to pan" },
  "hint.node": { zh: "节点: 点击空白处添加节点 (吸附 0.1m)", en: "Node: click empty area to add (snap 0.1m)" },
  "hint.element": { zh: "杆件: 依次点击两个节点创建单元, 点空白取消", en: "Member: click two nodes, click empty to cancel" },
  "hint.load": { zh: "荷载: 点击节点或杆件设置荷载", en: "Load: click a node or member to add load" },
  "hint.constraint": { zh: "约束: 点击节点, 在右侧勾选 ux/uy/rz", en: "Support: click a node, tick ux/uy/rz on the right" },
  "hint.erase": { zh: "删除: 点击荷载/约束/节点/杆件删除", en: "Erase: click load/support/node/member to delete" },

  // 右栏
  "llmAssistant": { zh: "LLM 建模助手", en: "LLM Assistant" },
  "settings": { zh: "设置", en: "Settings" },
  "solveResults": { zh: "求解结果", en: "Results" },
  "model": { zh: "模型", en: "Model" },
  "material": { zh: "材料", en: "Materials" },
  "section": { zh: "截面", en: "Sections" },
  "statusBar.ok": { zh: "✔ 模型有效", en: "✔ Model OK" },

  // 主题/语言/显示
  "theme": { zh: "主题", en: "Theme" },
  "theme.dark": { zh: "深色", en: "Dark" },
  "theme.light": { zh: "浅色", en: "Light" },
  "theme.ocean": { zh: "海洋", en: "Ocean" },
  "language": { zh: "语言", en: "Language" },
  "display": { zh: "显示", en: "Display" },
  "display.nodeLabels": { zh: "节点编号", en: "Node Labels" },
  "display.elementLabels": { zh: "杆件编号", en: "Member Labels" },
  "display.loads": { zh: "显示荷载", en: "Show Loads" },
  "display.constraints": { zh: "显示约束", en: "Show Supports" },

  // 荷载对话框
  "loadDialog.title": { zh: "添加荷载", en: "Add Load" },
  "loadDialog.atNode": { zh: "@节点 {n}", en: "@node {n}" },
  "loadDialog.type": { zh: "类型", en: "Type" },
  "loadDialog.direction": { zh: "方向", en: "Direction" },
  "loadDialog.dirX": { zh: "x (水平)", en: "x (horizontal)" },
  "loadDialog.dirY": { zh: "y (竖向)", en: "y (vertical)" },
  "loadDialog.dirR": { zh: "rz (转动)", en: "rz (rotation)" },
  "loadDialog.element": { zh: "作用单元", en: "Element" },
  "loadDialog.elemOption": { zh: "单元 #{id} (节点 {a}–{b})", en: "Elem #{id} (nodes {a}–{b})" },
  "loadDialog.t0": { zh: "下表面温变 T0 (℃)", en: "Bottom surface ΔT T0 (℃)" },
  "loadDialog.t1": { zh: "上表面温变 T1 (℃)", en: "Top surface ΔT T1 (℃)" },
  "loadDialog.value": { zh: "数值", en: "Value" },
  "loadDialog.valueDisp": { zh: "位移值 (m)", en: "Displacement (m)" },
  "loadDialog.valueUnit": { zh: "(N 或 N·m)", en: "(N or N·m)" },
  "loadDialog.valuePerM": { zh: "(N/m)", en: "(N/m)" },
  "loadDialog.position": { zh: "作用长度/位置 position (m, 从单元起点, 默认整跨)", en: "Length/position (m, from element start, default full span)" },
  "loadDialog.confirm": { zh: "确认", en: "OK" },
  "loadDialog.cancel": { zh: "取消", en: "Cancel" },
  "toast.needNodeLoad": { zh: "节点荷载需要先点击节点添加；当前是杆件，请选择单元荷载类型", en: "Nodal loads need a node: click a node first, or choose an element load type" },
  "toast.needElemLoad": { zh: "单元荷载需要先点击杆件添加；当前是节点，请选择节点荷载类型", en: "Element loads need a member: click a member first, or choose a nodal load type" },
  "toast.invalidElem": { zh: "请选择有效的单元", en: "Please select a valid element" },
  "toast.tempNeedsElem": { zh: "温度荷载需要选择作用单元", en: "Temperature load requires an element" },

  // 荷载类型
  "loadType.nodalForce": { zh: "节点集中力", en: "Nodal Force" },
  "loadType.lateralForce": { zh: "杆件横向集中力", en: "Lateral Point Load" },
  "loadType.lateralUniformPressure": { zh: "横向均布荷载", en: "Uniform Lateral Load" },
  "loadType.lateralLinearlyPressure": { zh: "横向线性分布荷载", en: "Linear Lateral Load" },
  "loadType.momentOnPoint": { zh: "节点弯矩", en: "Nodal Moment" },
  "loadType.axialForce": { zh: "杆件轴向集中力", en: "Axial Point Load" },
  "loadType.axialPressure": { zh: "杆件轴向均布荷载", en: "Uniform Axial Load" },
  "loadType.temperature": { zh: "温度荷载", en: "Temperature" },
  "loadType.supportMove": { zh: "支座位移", en: "Support Displacement" },

  // 画布标签/按钮
  "canvas.exportPng": { zh: "📷 导出图片", en: "📷 Export PNG" },
  "canvas.exportPngTip": { zh: "导出画布为 PNG 图片", en: "Export canvas as PNG" },
  "canvas.toast.exported": { zh: "已导出 PNG 图片", en: "PNG image exported" },
  "canvas.toast.exportFailed": { zh: "导出失败: 图片渲染错误", en: "Export failed: image render error" },
  "canvas.toast.exportErr": { zh: "导出失败: {msg}", en: "Export failed: {msg}" },
  "canvas.original": { zh: "原图", en: "Model" },
  "canvas.originalTip": { zh: "隐藏内力图", en: "Hide force diagrams" },
  "canvas.nDiagram": { zh: "轴力 N", en: "Axial N" },
  "canvas.vDiagram": { zh: "剪力 V", en: "Shear V" },
  "canvas.mDiagram": { zh: "弯矩 M", en: "Moment M" },
  "canvas.dDiagram": { zh: "挠度", en: "Deflection" },
  "canvas.nDiagramTip": { zh: "轴力图 N", en: "Axial force N" },
  "canvas.vDiagramTip": { zh: "剪力图 V", en: "Shear force V" },
  "canvas.mDiagramTip": { zh: "弯矩图 M", en: "Bending moment M" },
  "canvas.dDiagramTip": { zh: "挠度/变形图", en: "Deflection / deformed shape" },
  "canvas.udlLabel": { zh: "{v} 均布", en: "{v} UDL" },
  "canvas.axialLabel": { zh: "{v} 轴压", en: "{v} axial" },
  "canvas.linearLabel": { zh: "线性 {v} N/m", en: "Linear {v} N/m" },
  "canvas.tempLabel": { zh: "ΔT {a}~{b}℃", en: "ΔT {a}~{b}℃" },
  "canvas.dispUnit": { zh: "{v}mm", en: "{v}mm" },

  // 荷载/约束 Inspector
  "inspector.nodeTitle": { zh: "节点 #{id}", en: "Node #{id}" },
  "inspector.elemTitle": { zh: "杆件 #{id}", en: "Member #{id}" },
  "inspector.loadTitle": { zh: "荷载 #{id}", en: "Load #{id}" },
  "inspector.constraintTitle": { zh: "约束 @节点 {n}", en: "Support @node {n}" },
  "inspector.modelTitle": { zh: "模型", en: "Model" },
  "inspector.x": { zh: "x (m)", en: "x (m)" },
  "inspector.y": { zh: "y (m)", en: "y (m)" },
  "inspector.type": { zh: "类型", en: "Type" },
  "inspector.typeFrame": { zh: "frame (刚架)", en: "frame (rigid)" },
  "inspector.typeTruss": { zh: "truss (桁架)", en: "truss" },
  "inspector.section": { zh: "截面", en: "Section" },
  "inspector.material": { zh: "材料", en: "Material" },
  "inspector.node": { zh: "节点", en: "Node" },
  "inspector.element": { zh: "单元", en: "Element" },
  "inspector.direction": { zh: "方向", en: "Direction" },
  "inspector.value": { zh: "数值", en: "Value" },
  "inspector.dispValue": { zh: "位移值", en: "Displacement" },
  "inspector.position": { zh: "位置/长度", en: "Position/Length" },
  "inspector.deleteLoad": { zh: "删除荷载", en: "Delete Load" },
  "inspector.t0": { zh: "下表面 ΔT", en: "Bottom ΔT" },
  "inspector.t1": { zh: "上表面 ΔT", en: "Top ΔT" },
  "inspector.ux": { zh: "ux (水平)", en: "ux (horizontal)" },
  "inspector.uy": { zh: "uy (竖向)", en: "uy (vertical)" },
  "inspector.rz": { zh: "rz (转动)", en: "rz (rotation)" },
  "inspector.constraintHint": { zh: "勾选 = 该自由度被约束 (位移 0)", en: "Tick = DOF is constrained (displacement 0)" },
  "inspector.title": { zh: "标题", en: "Title" },
  "inspector.nodesElems": { zh: "节点 / 单元", en: "Nodes / Elements" },
  "inspector.E": { zh: "E (Pa)", en: "E (Pa)" },
  "inspector.mu": { zh: "μ", en: "μ" },
  "inspector.alpha": { zh: "α", en: "α" },
  "inspector.A": { zh: "A (m²)", en: "A (m²)" },
  "inspector.Iz": { zh: "Iz (m⁴)", en: "Iz (m⁴)" },
  "inspector.h": { zh: "h (m)", en: "h (m)" },

  // 结果面板
  "results.solveBtn": { zh: "求解", en: "Solve" },
  "results.solving": { zh: "求解中…", en: "Solving…" },
  "results.placeholder": { zh: "点击「求解」开始。示例: 2m 简支梁跨中 10kN。", en: "Click Solve to start. Example: 2m simply supported beam, 10kN mid-span." },
  "results.failed": { zh: "求解失败", en: "Solve failed" },
  "results.unknownErr": { zh: "未知错误", en: "Unknown error" },
  "results.statNodes": { zh: "节点", en: "Nodes" },
  "results.statElems": { zh: "单元", en: "Elements" },
  "results.statDof": { zh: "自由DOF", en: "Free DOF" },
  "results.deformScale": { zh: "变形倍率", en: "Deform scale" },
  "results.displacements": { zh: "位移 (m / rad)", en: "Displacements (m / rad)" },
  "results.endForces": { zh: "杆端力 (局部: N 轴力, V 剪力, M 弯矩)", en: "End Forces (local: N axial, V shear, M moment)" },
  "results.momentNote": { zh: "弯矩正负: M_i 逆时针为正; 画图时正弯矩在杆件上侧(受拉侧)。M 图跨节点连续。", en: "Moment sign: M_i positive counter-clockwise; positive moment drawn above the member (tension side). M diagram is continuous across nodes." },
  "results.reactions": { zh: "支座反力 (N / N·m)", en: "Reactions (N / N·m)" },
  "results.thNode": { zh: "节点", en: "Node" },
  "results.thElement": { zh: "单元", en: "Element" },
  "results.maxLabel": { zh: "max|{k}| (单元{e})", en: "max|{k}| (elem {e})" },

  // ChatPanel / LLM
  "chat.placeholder": { zh: "描述你的结构模型… (输入 /help 查看教程)", en: "Describe your structure… (/help for guide)" },
  "chat.send": { zh: "发送", en: "Send" },
  "chat.example": {
    zh: "🤖 我是结构力学建模 Agent，可以用自然语言帮你完成：\n\n① 建模：直接描述结构 → 我生成模型并应用\n　例：「建立 3x3 框架，每跨 1000mm，每层水平荷载 10kN」\n② 修改：在现有模型上调整\n　例：「把荷载改成 20kN」「加一层」\n③ 求解与解读：自动求解并汇报关键内力\n　例：「求一下支座反力」「最大弯矩在哪」\n\n📌 使用提示：\n· 描述尽量包含尺寸、跨度、层数、荷载大小与方向\n· 单位可用 mm / kN（我会自动换算）\n· 输入 /help 可查看命令与内置教程工程\n· 输入 /new 新建空白项目\n\n试试直接说：「建一个 2m 简支梁，跨中 10kN 竖向荷载」",
    en: "🤖 I'm the structural modeling Agent. I can help you in natural language:\n\n① Modeling: describe a structure → I generate & apply it\n　e.g. \"Build a 3x3 frame, 1000mm per bay, 10kN lateral load per floor\"\n② Editing: adjust the current model\n　e.g. \"Change the load to 20kN\" \"Add one more story\"\n③ Solving & explanation: auto-solve and report key forces\n　e.g. \"What are the support reactions?\" \"Where is the max moment?\"\n\n📌 Tips:\n· Include dimensions, spans, stories, load values & directions\n· Units: mm / kN are fine (I'll convert)\n· Type /help for commands & built-in tutorial projects\n· Type /new to create a blank project\n\nTry: \"Build a 2m simply supported beam with a 10kN mid-span load\""
  },
  "chat.needsConfig": { zh: "请先点击「设置」配置 base_url / api_key / model。", en: "Please configure base_url / api_key / model in Settings first." },
  "chat.sameModel": { zh: "⚠️ 模型似乎与当前画布完全相同 —— LLM 可能没有按你的要求生成新模型，请重新描述需求，或先新建一个空白项目再提问。", en: "⚠️ The model seems identical to the current canvas — the LLM may not have generated a new model per your request. Please re-describe, or start a new project first." },
  "chat.requestFailed": { zh: "请求失败: {msg}", en: "Request failed: {msg}" },
  "chat.applied": { zh: "已应用模型（{n} 节点 / {e} 单元），正在求解…", en: "Model applied ({n} nodes / {e} elements), solving…" },
  "chat.configSaved": { zh: "LLM 配置已保存", en: "LLM config saved" },
  "chat.configSaveFailed": { zh: "保存配置失败: {msg}", en: "Save config failed: {msg}" },
  "chat.saveConfig": { zh: "保存配置", en: "Save Config" },
  "chat.generated": { zh: "已生成模型: {n} 节点 / {e} 单元", en: "Model generated: {n} nodes / {e} elements" },
  "chat.apply": { zh: "应用到画布 + 求解", en: "Apply + Solve" },
  "chat.settings": { zh: "设置", en: "Settings" },

  // store toast
  "toast.noUndo": { zh: "没有可撤销的操作", en: "Nothing to undo" },
  "toast.noRedo": { zh: "没有可重做的操作", en: "Nothing to redo" },
  "toast.solved": { zh: "求解完成", en: "Solved" },
  "toast.solveFailed": { zh: "求解失败: {msg}", en: "Solve failed: {msg}" },
  "toast.backendFailed": { zh: "调用后端失败: {msg}", en: "Backend error: {msg}" },

  // model.js 校验消息
  "validate.empty": { zh: "模型为空", en: "Model is empty" },
  "validate.minNodes": { zh: "至少需要 2 个节点", en: "At least 2 nodes required" },
  "validate.minElems": { zh: "至少需要 1 个单元", en: "At least 1 element required" },
  "validate.minCons": { zh: "至少需要 1 个约束", en: "At least 1 constraint required" },
  "validate.badElem": { zh: "单元 {id} 引用了不存在的节点", en: "Element {id} references a missing node" },
  "validate.badLoad": { zh: "荷载引用了不存在的节点", en: "Load references a missing node" },
  "validate.badCons": { zh: "约束引用了不存在的节点", en: "Constraint references a missing node" },
  "validate.isolatedNode": { zh: "节点 {id} 未连接任何单元（孤立节点），请删除或连接它", en: "Node {id} is not connected to any element (isolated node). Delete or connect it" },
  "validate.llmDefault": { zh: "LLM 生成模型", en: "LLM Generated Model" },

  // ipc.js
  "ipc.notInTauri": { zh: "当前不在 Tauri 环境中(请通过桌面应用打开)。", en: "Not running in Tauri (please open via the desktop app)." },

  // /help 与教程工程
  "help.intro": { zh: "🛠 可用命令：\n• /help — 显示本帮助\n• /new — 新建空白项目\n\n点击下方教程工程卡片可一键加载示例模型（含求解）。", en: "🛠 Commands:\n• /help — show this help\n• /new — create a blank project\n\nClick a tutorial project below to load an example model (solve included)." },
  "help.newDone": { zh: "已新建空白项目。", en: "Blank project created." },
  "help.tutorialTitle": { zh: "📚 教程工程", en: "📚 Tutorial Projects" },
  "help.loadTip": { zh: "加载此示例到画布", en: "Load this example to canvas" },
  "help.loaded": { zh: "已加载教程工程: {name}", en: "Tutorial loaded: {name}" },
  "help.loadFailed": { zh: "加载教程工程失败: {msg}", en: "Failed to load tutorial: {msg}" },
};

let currentLang = "zh";

export function setLang(lang) {
  currentLang = lang === "en" ? "en" : "zh";
  document.documentElement.lang = currentLang;
}

export function getLang() {
  return currentLang;
}

export function t(key, vars) {
  const entry = dict[key];
  if (!entry) return key;
  let s = entry[currentLang] || entry.zh;
  if (vars) {
    for (const [k, v] of Object.entries(vars)) {
      s = s.split("{" + k + "}").join(String(v));
    }
  }
  return s;
}

export const LANG_OPTIONS = [
  { id: "zh", label: "中文" },
  { id: "en", label: "English" },
];
