// 模型 → SVG 元素几何计算（GUI_REDESIGN_PLAN §4.6）
// 返回描述对象数组, 由 CanvasView 渲染成 JSX/SVG 元素。

import { worldToScreen } from "./transform.js";

/** 节点: circle + 编号 text 的屏幕坐标 */
export function computeNodes(model, view, selection) {
  return (model.nodes || []).map((n) => {
    const p = worldToScreen(n.x, n.y, view);
    return {
      id: n.id,
      cx: p.x,
      cy: p.y,
      type: n.type,
      selected: selection && selection.type === "node" && selection.id === n.id,
    };
  });
}

/** 杆件: line 两端屏幕坐标; truss 虚线, frame 实线 */
export function computeElements(model, view, selection) {
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  return (model.elements || [])
    .map((e) => {
      const a = byId.get(e.nodeI);
      const b = byId.get(e.nodeJ);
      if (!a || !b) return null;
      const pa = worldToScreen(a.x, a.y, view);
      const pb = worldToScreen(b.x, b.y, view);
      return {
        id: e.id,
        x1: pa.x,
        y1: pa.y,
        x2: pb.x,
        y2: pb.y,
        type: e.type,
        selected: selection && selection.type === "element" && selection.id === e.id,
      };
    })
    .filter(Boolean);
}

/**
 * 荷载: 从节点出发的箭头 (line + arrowhead polygon), 旁标数值。
 * direction x/y → 直线箭头; rz → 绕节点的小圆弧箭头。
 */
export function computeLoads(model, view, selection) {
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  return (model.loads || [])
    .map((ld) => {
      const n = byId.get(ld.node);
      if (!n) return null;
      const p = worldToScreen(n.x, n.y, view);
      const base = {
        id: ld.id,
        node: ld.node,
        direction: ld.direction,
        value: ld.value,
        x: p.x,
        y: p.y,
        selected: selection && selection.type === "load" && selection.id === ld.id,
      };

      if (ld.direction === "rz") {
        // 小圆弧箭头: 半径 16, 起于 0°, 按 value 符号决定扫掠方向
        const R = 16;
        const aEnd = ((ld.value >= 0 ? 235 : 125) * Math.PI) / 180;
        const sx0 = p.x + R;
        const sy0 = p.y;
        const ex = p.x + R * Math.cos(aEnd);
        const ey = p.y - R * Math.sin(aEnd);
        const sweep = ld.value >= 0 ? 1 : 0;
        const d = `M ${sx0},${sy0} A ${R},${R} 0 0 ${sweep} ${ex},${ey}`;
        const ux = -Math.sin(aEnd);
        const uy = -Math.cos(aEnd); // 弧端点切线方向
        const nrm = Math.hypot(ux, uy) || 1;
        const px = -uy / nrm;
        const py = ux / nrm;
        const hx = ex - (ux / nrm) * 11;
        const hy = ey - (uy / nrm) * 11;
        return {
          ...base,
          kind: "rz",
          d,
          head: { x: ex, y: ey },
          arrowPoints: `${ex},${ey} ${hx - px * 5},${hy - py * 5} ${hx + px * 5},${hy + py * 5}`,
        };
      }

      // 直线箭头: 力的方向 = value 符号; y 方向在屏幕上需翻转
      let ux = 0;
      let uy = 0;
      if (ld.direction === "x") ux = Math.sign(ld.value);
      else if (ld.direction === "y") uy = -Math.sign(ld.value);
      const L = 28;
      const tx = p.x + ux * L;
      const ty = p.y + uy * L;
      const nrm = Math.hypot(ux, uy) || 1;
      const dx = ux / nrm;
      const dy = uy / nrm;
      const px = -dy;
      const py = dx;
      const hx = tx - dx * 11;
      const hy = ty - dy * 11;
      return {
        ...base,
        kind: "arrow",
        tx,
        ty,
        arrowPoints: `${tx},${ty} ${hx - px * 5},${hy - py * 5} ${hx + px * 5},${hy + py * 5}`,
        labelX: (p.x + tx) / 2 + 7,
        labelY: (p.y + ty) / 2 - 5,
      };
    })
    .filter(Boolean);
}

/** 约束: 节点下方三角符号; ux+uy=固定(实心), 单平移=滚轴(空心+圆), rz=虚线圆弧 */
export function computeConstraints(model, view, selection) {
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  return (model.constraints || [])
    .map((c) => {
      const n = byId.get(c.node);
      if (!n) return null;
      const p = worldToScreen(n.x, n.y, view);
      const h = 11;
      const half = 9;
      const y = p.y + 9;
      return {
        node: c.node,
        dofs: c.dofs,
        x: p.x,
        y,
        h,
        half,
        pts: `${p.x},${y} ${p.x - half},${y + h} ${p.x + half},${y + h}`,
        fixed: c.dofs.includes("ux") && c.dofs.includes("uy"),
        roller: c.dofs.length === 1,
        hasRz: c.dofs.includes("rz"),
        selected: selection && selection.type === "constraint" && selection.id === c.node,
      };
    })
    .filter(Boolean);
}

/** 变形图: newPos = pos + u*deformScale, 红色虚线叠加 */
export function computeDeformed(model, results, deformScale, view) {
  if (!results || !Array.isArray(results.displacements)) return [];
  const disp = new Map(results.displacements.map((d) => [d.node, d]));
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  return (model.elements || [])
    .map((e) => {
      const a = byId.get(e.nodeI);
      const b = byId.get(e.nodeJ);
      if (!a || !b) return null;
      const da = disp.get(e.nodeI);
      const db = disp.get(e.nodeJ);
      const ax = a.x + (da ? (da.ux || 0) * deformScale : 0);
      const ay = a.y + (da ? (da.uy || 0) * deformScale : 0);
      const bx = b.x + (db ? (db.ux || 0) * deformScale : 0);
      const by = b.y + (db ? (db.uy || 0) * deformScale : 0);
      const pa = worldToScreen(ax, ay, view);
      const pb = worldToScreen(bx, by, view);
      return { x1: pa.x, y1: pa.y, x2: pb.x, y2: pb.y };
    })
    .filter(Boolean);
}

/** 自动估算变形倍率: 0.15*span/max|u|; max|u|=0 时返回 0 */
export function estimateDeformScale(model, results) {
  if (!results || !Array.isArray(results.displacements) || !results.displacements.length) return 0;
  let maxU = 0;
  for (const d of results.displacements) {
    maxU = Math.max(maxU, Math.abs(d.ux || 0), Math.abs(d.uy || 0));
  }
  if (maxU < 1e-12) return 0;
  let minX = Infinity,
    maxX = -Infinity;
  for (const n of model.nodes || []) {
    minX = Math.min(minX, n.x);
    maxX = Math.max(maxX, n.x);
  }
  const span = Math.max(maxX - minX, 1);
  return (0.15 * span) / maxU;
}

/**
 * 内力图: 在杆件上叠加 N/V/M 图 (局部坐标端力 → 屏幕偏移折线)
 * kind: "N" | "V" | "M"
 * 每根杆: 局部端值 fI, fJ → 垂直杆轴偏移(米) → 屏幕折线
 * N 图沿杆轴方向偏移(画在杆两侧), V/M 垂直杆轴偏移(画在杆一侧, 正负分侧)
 */
export function computeForceDiagram(model, results, kind, view) {
  if (!results || !Array.isArray(results.endForces)) return [];
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  const forceById = new Map(results.endForces.map((f) => [f.element, f]));
  const pickSide = (f) => (kind === "N" ? "N" : kind === "V" ? "V" : "M");

  // 全局最大 |值| 用于归一化
  let maxVal = 1e-9;
  const side = pickSide(null);
  for (const f of results.endForces) {
    const ni = f.nodeI || {};
    const nj = f.nodeJ || {};
    maxVal = Math.max(maxVal, Math.abs(ni[side] || 0), Math.abs(nj[side] || 0));
  }

  const out = [];
  for (const el of model.elements || []) {
    const a = byId.get(el.nodeI);
    const b = byId.get(el.nodeJ);
    const f = forceById.get(el.id);
    if (!a || !b || !f) continue;

    // 杆轴方向 (世界坐标)
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.hypot(dx, dy) || 1e-9;
    const ux = dx / len;
    const uy = dy / len;
    // 垂直方向: 屏幕 y 向下, 世界 y 向上 → 取 (uy, -ux), 使正 V/M 画在视觉上方
    const nx = uy;
    const ny = -ux;

    // 偏移量(世界单位, 米): 最大偏移 = 0.25 * 杆长
    const amp = 0.25 * len;
    const k = amp / maxVal;
    const ni = f.nodeI || {};
    const nj = f.nodeJ || {};
    const fI = ni[side] || 0;
    const fJ = nj[side] || 0;
    const oI = fI * k;
    const oJ = fJ * k;

    // 折线采样: 直线段两端的屏幕点
    const pI = worldToScreen(a.x + nx * oI, a.y + ny * oI, view);
    const pJ = worldToScreen(b.x + nx * oJ, b.y + ny * oJ, view);
    out.push({
      key: `fg-${kind}-${el.id}`,
      d: `M ${pI.x},${pI.y} L ${pJ.x},${pJ.y}`,
    });
  }
  return out;
}
