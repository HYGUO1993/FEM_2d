import { create } from "zustand";
import { defaultModel, newEmptyModel, nextId, validateModel } from "./model.js";
import { estimateDeformScale } from "./canvas/render.js";
import * as ipc from "./ipc.js";

let toastTimer = null;

// UI 布局偏好 (localStorage 持久化)
const savedUi = (() => {
  try {
    return JSON.parse(localStorage.getItem("femlab_ui") || "null");
  } catch {
    return null;
  }
})();

const initialState = {
  projects: [],
  currentProject: "",
  model: defaultModel(),
  tool: "select", // select | node | element | load | constraint | erase
  selection: null, // {type:'node'|'element'|'load'|'constraint', id}
  results: null,
  solved: false,
  pendingNodeA: null,
  llmConfig: { base_url: "", api_key: "", model: "" },
  chatMessages: [],
  deformScale: 30,
  view: { scale: 50, panX: 0, panY: 0 }, // {scale, panX, panY} 非 null, 避免首帧渲染崩溃
  viewBox: null, // {w, h}
  toast: null, // {msg, isError}
  loadDialog: null, // {node}
  pendingLlmModel: null,
  solveTime: "",
  solving: false,
  // —— 布局偏好 ——
  sidebarWidth: savedUi?.sidebarWidth ?? 240, // 左栏宽 px
  rightWidth: savedUi?.rightWidth ?? 320,     // 右栏宽 px
  llmPosition: savedUi?.llmPosition ?? "right", // "right" | "bottom"
};

export const useStore = create((set, get) => ({
  ...initialState,

  // —— 项目 ——
  setProjects: (projects) => set({ projects }),
  setCurrentProject: (name) => set({ currentProject: name }),

  // —— 工具 / 选择 ——
  setTool: (tool) => set({ tool, pendingNodeA: null }),
  select: (sel) => set({ selection: sel }),

  // 整体替换模型（加载项目 / LLM 应用）
  setModel: (model) =>
    set({ model, selection: null, results: null, solved: false, pendingNodeA: null }),

  // —— 节点 ——
  addNode: (x, y) =>
    set((s) => {
      const id = nextId(s.model.nodes);
      const nodes = [...s.model.nodes, { id, type: "frame", x, y }];
      return {
        model: { ...s.model, nodes },
        results: null,
        solved: false,
        selection: { type: "node", id },
      };
    }),

  moveNode: (id, x, y) =>
    set((s) => ({
      model: { ...s.model, nodes: s.model.nodes.map((n) => (n.id === id ? { ...n, x, y } : n)) },
      results: null,
      solved: false,
    })),

  updateNode: (id, patch) =>
    set((s) => ({
      model: { ...s.model, nodes: s.model.nodes.map((n) => (n.id === id ? { ...n, ...patch } : n)) },
      results: null,
      solved: false,
    })),

  // 级联删除: 引用该节点的 elements/loads/constraints 一并删除
  deleteNode: (id) =>
    set((s) => {
      const m = s.model;
      return {
        model: {
          ...m,
          nodes: m.nodes.filter((n) => n.id !== id),
          elements: m.elements.filter((e) => e.nodeI !== id && e.nodeJ !== id),
          loads: m.loads.filter((l) => l.node !== id),
          constraints: m.constraints.filter((c) => c.node !== id),
        },
        results: null,
        solved: false,
        selection: null,
        pendingNodeA: s.pendingNodeA === id ? null : s.pendingNodeA,
      };
    }),

  // —— 单元 ——
  addElement: (nodeI, nodeJ) =>
    set((s) => {
      if (nodeI == null || nodeJ == null || nodeI === nodeJ) return { pendingNodeA: null };
      const dup = s.model.elements.some(
        (e) => (e.nodeI === nodeI && e.nodeJ === nodeJ) || (e.nodeI === nodeJ && e.nodeJ === nodeI)
      );
      if (dup) return { pendingNodeA: null };
      const id = nextId(s.model.elements);
      const elements = [
        ...s.model.elements,
        { id, type: "frame", nodeI, nodeJ, section: 0, material: 0 },
      ];
      return {
        model: { ...s.model, elements },
        results: null,
        solved: false,
        pendingNodeA: null,
      };
    }),

  updateElement: (id, patch) =>
    set((s) => ({
      model: { ...s.model, elements: s.model.elements.map((e) => (e.id === id ? { ...e, ...patch } : e)) },
      results: null,
      solved: false,
    })),

  deleteElement: (id) =>
    set((s) => ({
      model: { ...s.model, elements: s.model.elements.filter((e) => e.id !== id) },
      results: null,
      solved: false,
      selection: null,
    })),

  // —— 荷载 ——
  addLoad: (node, direction, value) =>
    set((s) => {
      const id = nextId(s.model.loads);
      const loads = [...s.model.loads, { type: "nodalForce", direction, value, node }];
      return {
        model: { ...s.model, loads },
        results: null,
        solved: false,
        loadDialog: null,
      };
    }),

  deleteLoad: (id) =>
    set((s) => ({
      model: { ...s.model, loads: s.model.loads.filter((l) => l.id !== id) },
      results: null,
      solved: false,
      selection: null,
    })),

  // —— 约束 ——
  toggleConstraint: (node, dof) =>
    set((s) => {
      const m = s.model;
      const existing = m.constraints.find((c) => c.node === node);
      let constraints;
      if (existing) {
        const has = existing.dofs.includes(dof);
        const dofs = has ? existing.dofs.filter((d) => d !== dof) : [...existing.dofs, dof];
        constraints = dofs.length
          ? m.constraints.map((c) => (c.node === node ? { ...c, dofs } : c))
          : m.constraints.filter((c) => c.node !== node);
      } else {
        constraints = [...m.constraints, { node, dofs: [dof] }];
      }
      return {
        model: { ...m, constraints },
        results: null,
        solved: false,
      };
    }),

  removeConstraint: (node) =>
    set((s) => ({
      model: { ...s.model, constraints: s.model.constraints.filter((c) => c.node !== node) },
      results: null,
      solved: false,
    })),

  // —— 材料 / 截面 ——
  updateMaterial: (id, patch) =>
    set((s) => ({
      model: {
        ...s.model,
        materials: s.model.materials.map((mt) => (mt.id === id ? { ...mt, ...patch } : mt)),
      },
      results: null,
      solved: false,
    })),

  updateSection: (id, patch) =>
    set((s) => ({
      model: {
        ...s.model,
        sections: s.model.sections.map((sc) => (sc.id === id ? { ...sc, ...patch } : sc)),
      },
      results: null,
      solved: false,
    })),

  // —— 视图 ——
  setView: (view) => set({ view }),
  setViewBox: (w, h) => set({ viewBox: { w, h } }),
  resetView: () => {
    const { viewBox } = get();
    if (viewBox) get().fitView(viewBox.w, viewBox.h);
  },

  // —— 布局偏好 ——
  setSidebarWidth: (w) => {
    set({ sidebarWidth: w });
    try {
      localStorage.setItem(
        "femlab_ui",
        JSON.stringify({ sidebarWidth: w, rightWidth: get().rightWidth, llmPosition: get().llmPosition })
      );
    } catch {}
  },
  setRightWidth: (w) => {
    set({ rightWidth: w });
    try {
      localStorage.setItem(
        "femlab_ui",
        JSON.stringify({ sidebarWidth: get().sidebarWidth, rightWidth: w, llmPosition: get().llmPosition })
      );
    } catch {}
  },
  setLlmPosition: (pos) => {
    set({ llmPosition: pos });
    try {
      localStorage.setItem(
        "femlab_ui",
        JSON.stringify({ sidebarWidth: get().sidebarWidth, rightWidth: get().rightWidth, llmPosition: pos })
      );
    } catch {}
  },

  fitView: (w, h) => {
    const m = get().model;
    if (!m.nodes.length) {
      set({ view: { scale: 50, panX: w / 2, panY: h / 2 } });
      return;
    }
    let minX = Infinity,
      maxX = -Infinity,
      minY = Infinity,
      maxY = -Infinity;
    for (const n of m.nodes) {
      minX = Math.min(minX, n.x);
      maxX = Math.max(maxX, n.x);
      minY = Math.min(minY, n.y);
      maxY = Math.max(maxY, n.y);
    }
    const spanX = Math.max(maxX - minX, 0.5);
    const spanY = Math.max(maxY - minY, 0.5);
    const scale = Math.min((w * 0.6) / spanX, (h * 0.6) / spanY);
    const cx = (minX + maxX) / 2;
    const cy = (minY + maxY) / 2;
    set({ view: { scale, panX: w / 2 - cx * scale, panY: h / 2 + cy * scale } });
  },

  // —— 结果 ——
  setResults: (result, solveTimeStr) => {
    const m = get().model;
    let def = 30;
    if (result && result.status === "ok") {
      const est = estimateDeformScale(m, result);
      if (est > 0) def = est;
    }
    set({
      results: result,
      solved: !!(result && result.status === "ok"),
      deformScale: def,
      solveTime: solveTimeStr || "",
    });
  },

  async solve() {
    const s = get();
    const err = validateModel(s.model);
    if (err) {
      set({ toast: { msg: err, isError: true } });
      return null;
    }
    set({ solving: true, solved: false });
    const t0 = performance.now();
    try {
      const resultStr = await ipc.solve(s.model);
      const result = JSON.parse(resultStr);
      const secs = ((performance.now() - t0) / 1000).toFixed(2);
      get().setResults(result, secs);
      set({
        toast: {
          msg: result.status === "ok" ? "求解完成" : "求解失败: " + (result.message || ""),
          isError: result.status !== "ok",
        },
      });
      return result;
    } catch (e) {
      const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
      set({ toast: { msg: "调用后端失败: " + msg, isError: true } });
      return null;
    } finally {
      set({ solving: false });
    }
  },

  // —— LLM ——
  setLlmConfig: (config) => set({ llmConfig: config }),
  pushChat: (role, content) => set((s) => ({ chatMessages: [...s.chatMessages, { role, content }] })),
  setPendingLlmModel: (m) => set({ pendingLlmModel: m }),

  // —— UI ——
  setDeformScale: (v) => set({ deformScale: v }),
  setToast: (msg, isError = false) => {
    clearTimeout(toastTimer);
    set({ toast: { msg, isError } });
    toastTimer = setTimeout(() => set({ toast: null }), 2600);
  },
  clearToast: () => set({ toast: null }),
  setLoadDialog: (d) => set({ loadDialog: d }),
  setPendingNodeA: (id) => set({ pendingNodeA: id }),
}));
