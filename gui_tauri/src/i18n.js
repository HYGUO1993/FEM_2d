// 轻量 i18n: zh / en 字典 + t() 函数
// 用法: import { t, setLang } from "./i18n.js"

const dict = {
  // 顶栏
  "llm.toBottom": { zh: "LLM 移到底部", en: "LLM to Bottom" },
  "llm.toRight": { zh: "LLM 移到右侧", en: "LLM to Right" },
  "undo": { zh: "撤销", en: "Undo" },
  "redo": { zh: "重做", en: "Redo" },
  "resetView": { zh: "重置视图", en: "Reset View" },
  "solve": { zh: "求解", en: "Solve" },
  "solving": { zh: "求解中…", en: "Solving…" },
  "undoTip": { zh: "撤销 (Ctrl+Z)", en: "Undo (Ctrl+Z)" },
  "redoTip": { zh: "重做 (Ctrl+Y)", en: "Redo (Ctrl+Y)" },
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
  "statusBar.err": { zh: "✘ ", en: "✘ " },
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
};

let currentLang = "zh";

export function setLang(lang) {
  currentLang = lang === "en" ? "en" : "zh";
  document.documentElement.lang = currentLang;
}

export function getLang() {
  return currentLang;
}

export function t(key) {
  const entry = dict[key];
  if (!entry) return key;
  return entry[currentLang] || entry.zh;
}

export const LANG_OPTIONS = [
  { id: "zh", label: "中文" },
  { id: "en", label: "English" },
];
