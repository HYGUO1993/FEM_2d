// FEM 模型 schema（见 GUI_REDESIGN_PLAN §2）与工厂函数
// 字段名与 femcli 严格一致: nodes/constraints/elements/materials/sections/loads,
// nodeI/nodeJ/dofs/nodalForce 等, 不得改动。

// 荷载类型 (与 femcli.cpp LoadTypeFromStr 一致)
export const LOAD_TYPES = {
  nodalForce: "节点集中力",
  lateralForce: "杆件横向集中力",
  lateralUniformPressure: "横向均布荷载",
  lateralLinearlyPressure: "横向线性分布荷载",
  momentOnPoint: "节点弯矩",
  axialForce: "杆件轴向集中力",
  axialPressure: "杆件轴向均布荷载",
  temperature: "温度荷载",
  supportMove: "支座位移",
};

export function defaultModel() {
  return {
    schemaVersion: "1.0",
    title: "2m 简支梁跨中集中力",
    solver: "builtin",
    nodes: [
      { id: 0, type: "frame", x: 0.0, y: 0.0 },
      { id: 1, type: "frame", x: 1.0, y: 0.0 },
      { id: 2, type: "frame", x: 2.0, y: 0.0 },
    ],
    constraints: [
      { node: 0, dofs: ["ux", "uy"] },
      { node: 2, dofs: ["uy"] },
    ],
    elements: [
      { id: 0, type: "frame", nodeI: 0, nodeJ: 1, section: 0, material: 0 },
      { id: 1, type: "frame", nodeI: 1, nodeJ: 2, section: 0, material: 0 },
    ],
    materials: [{ id: 0, E: 210000000000.0, mu: 0.3, alpha: 0.0 }],
    sections: [{ id: 0, A: 0.01, Iz: 1e-5, height: 0.1 }],
    loads: [{ type: "nodalForce", direction: "y", value: -10000.0, node: 1 }],
  };
}

export function defaultMaterial(id = 0) {
  return { id, E: 210000000000.0, mu: 0.3, alpha: 0.0 };
}

export function defaultSection(id = 0) {
  return { id, A: 0.01, Iz: 1e-5, height: 0.1 };
}

export function newEmptyModel(title) {
  return {
    schemaVersion: "1.0",
    title: title || "未命名项目",
    solver: "builtin",
    nodes: [],
    constraints: [],
    elements: [],
    materials: [defaultMaterial(0)],
    sections: [defaultSection(0)],
    loads: [],
  };
}

/** id 自增: 当前最大 + 1（空数组返回 0） */
export function nextId(arr) {
  return arr.reduce((m, o) => Math.max(m, o.id ?? -1), -1) + 1;
}

/** 节点自由度数: truss=2, frame=3 */
export function nodeDofCount(node) {
  return node.type === "truss" ? 2 : 3;
}

/** 前端校验: 通过返回 null, 否则返回错误字符串 */
export function validateModel(m) {
  if (!m) return "模型为空";
  if (!Array.isArray(m.nodes) || m.nodes.length < 2) return "至少需要 2 个节点";
  if (!Array.isArray(m.elements) || m.elements.length < 1) return "至少需要 1 个单元";
  if (!Array.isArray(m.constraints) || m.constraints.length < 1) return "至少需要 1 个约束";
  const ids = new Set(m.nodes.map((n) => n.id));
  for (const e of m.elements) {
    if (!ids.has(e.nodeI) || !ids.has(e.nodeJ)) return `单元 ${e.id} 引用了不存在的节点`;
  }
  for (const l of m.loads || []) {
    if (!ids.has(l.node)) return "荷载引用了不存在的节点";
  }
  for (const c of m.constraints) {
    if (!ids.has(c.node)) return "约束引用了不存在的节点";
  }
  return null;
}

export function modelSummary(m) {
  return {
    nodes: (m.nodes || []).length,
    elements: (m.elements || []).length,
    dofs: (m.nodes || []).reduce((acc, n) => acc + nodeDofCount(n), 0),
  };
}

/** 归一化 LLM 返回的模型: 补齐缺失数组与默认材料/截面 */
export function normalizeModel(m) {
  const out = {
    schemaVersion: "1.0",
    title: m.title || "LLM 生成模型",
    solver: "builtin",
    nodes: Array.isArray(m.nodes) ? m.nodes : [],
    constraints: Array.isArray(m.constraints) ? m.constraints : [],
    elements: Array.isArray(m.elements) ? m.elements : [],
    materials: Array.isArray(m.materials) && m.materials.length ? m.materials : [defaultMaterial(0)],
    sections: Array.isArray(m.sections) && m.sections.length ? m.sections : [defaultSection(0)],
    loads: Array.isArray(m.loads) ? m.loads : [],
  };
  // 保证节点/单元/材料/截面引用到的 id 都存在
  out.elements = out.elements.map((e) => ({
    ...e,
    section: e.section ?? 0,
    material: e.material ?? 0,
  }));
  return out;
}
