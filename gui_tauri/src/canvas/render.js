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
 * 荷载: 支持节点集中力/弯矩 + 杆件集中力/均布/线性分布/温度
 * 返回描述对象数组, CanvasView 按 kind 渲染
 */
export function computeLoads(model, view, selection) {
  const nodeById = new Map((model.nodes || []).map((n) => [n.id, n]));
  const elemById = new Map((model.elements || []).map((e) => [e.id, e]));

  // 世界坐标 → 屏幕
  const P = (x, y) => worldToScreen(x, y, view);

  return (model.loads || [])
    .map((ld) => {
      const sel = selection && selection.type === "load" && selection.id === ld.id;
      const base = {
        id: ld.id,
        node: ld.node,
        element: ld.element,
        type: ld.type,
        direction: ld.direction,
        value: ld.value,
        selected: sel,
      };

      // —— 温度荷载: 杆上红色温度标注 ——
      if (ld.type === "temperature") {
        const e = elemById.get(ld.element);
        if (!e) return null;
        const a = nodeById.get(e.nodeI);
        const b = nodeById.get(e.nodeJ);
        if (!a || !b) return null;
        const pa = P(a.x, a.y);
        const pb = P(b.x, b.y);
        const mx = (pa.x + pb.x) / 2;
        const my = (pa.y + pb.y) / 2;
        return {
          ...base,
          kind: "temp",
          x: mx,
          y: my - 12,
          label: `ΔT ${ld.T0 ?? 0}~${ld.T1 ?? 0}℃`,
        };
      }

      // —— 杆件单元荷载 ——
      if (ld.type === "lateralForce" || ld.type === "lateralUniformPressure" ||
          ld.type === "lateralLinearlyPressure" || ld.type === "axialForce" ||
          ld.type === "axialPressure") {
        const e = elemById.get(ld.element);
        if (!e) return null;
        const a = nodeById.get(e.nodeI);
        const b = nodeById.get(e.nodeJ);
        if (!a || !b) return null;
        const pa = P(a.x, a.y);
        const pb = P(b.x, b.y);
        // 杆轴单位向量 (屏幕)
        const len = Math.hypot(pb.x - pa.x, pb.y - pa.y) || 1;
        const ux = (pb.x - pa.x) / len;
        const uy = (pb.y - pa.y) / len;
        // 垂直方向 (局部 y): 屏幕向下, 取 (uy, -ux)
        const vx = uy;
        const vy = -ux;

        // 力的方向: 横向(direction=x/y) → 垂直杆轴; 轴向 → 沿杆轴
        let fx, fy; // 力的单位方向(屏幕)
        if (ld.type === "axialForce" || ld.type === "axialPressure") {
          const s = ld.direction === "x" ? 1 : -1;
          fx = ux * s * Math.sign(ld.value || 1);
          fy = uy * s * Math.sign(ld.value || 1);
        } else {
          // 横向: 屏幕垂直方向, value 符号决定正负
          const s = Math.sign(ld.value || 1);
          fx = vx * s;
          fy = vy * s;
        }

        // 分布荷载: 沿线画 n 个箭头
        if (ld.type === "lateralUniformPressure" || ld.type === "axialPressure") {
          const n = Math.max(3, Math.floor(len / 22));
          const arrows = [];
          for (let i = 1; i <= n; i++) {
            const t = i / (n + 1);
            const sx = pa.x + (pb.x - pa.x) * t;
            const sy = pa.y + (pb.y - pa.y) * t;
            arrows.push({ sx, sy, ...arrowGeom(sx, sy, fx, fy, 16) });
          }
          return {
            ...base,
            kind: "udl",
            arrows,
            label: `${fmtLoad(ld.value)} ${ld.type === "axialPressure" ? "轴压" : "均布"}`,
            labelX: (pa.x + pb.x) / 2 + 6,
            labelY: (pa.y + pb.y) / 2 - 14,
          };
        }
        // 线性分布: 箭头长度从 0 → max
        if (ld.type === "lateralLinearlyPressure") {
          const n = Math.max(3, Math.floor(len / 22));
          const arrows = [];
          for (let i = 1; i <= n; i++) {
            const t = i / (n + 1);
            const sx = pa.x + (pb.x - pa.x) * t;
            const sy = pa.y + (pb.y - pa.y) * t;
            const mag = (t * 22).toFixed(1);
            arrows.push({ sx, sy, ...arrowGeom(sx, sy, fx, fy, Number(mag)) });
          }
          return {
            ...base,
            kind: "udl",
            arrows,
            label: `线性 ${fmtLoad(ld.value)} N/m`,
            labelX: (pa.x + pb.x) / 2 + 6,
            labelY: (pa.y + pb.y) / 2 - 14,
          };
        }
        // 杆上集中力
        const pos = ld.position || 0;
        const t = pos / (Math.hypot(b.x - a.x, b.y - a.y) || 1);
        const tClamped = Math.max(0, Math.min(1, t));
        const sx = pa.x + (pb.x - pa.x) * tClamped;
        const sy = pa.y + (pb.y - pa.y) * tClamped;
        const g = arrowGeom(sx, sy, fx, fy, 24);
        return {
          ...base,
          kind: "elemArrow",
          sx,
          sy,
          ...g,
          labelX: sx + fx * 30 + 6,
          labelY: sy + fy * 30 - 4,
        };
      }

      // —— 节点荷载 ——
      const n = nodeById.get(ld.node);
      if (!n) return null;
      const p = P(n.x, n.y);

      if (ld.direction === "rz" || ld.type === "momentOnPoint") {
        const R = 16;
        const aEnd = ((ld.value >= 0 ? 235 : 125) * Math.PI) / 180;
        const sx0 = p.x + R;
        const sy0 = p.y;
        const ex = p.x + R * Math.cos(aEnd);
        const ey = p.y - R * Math.sin(aEnd);
        const sweep = ld.value >= 0 ? 1 : 0;
        const d = `M ${sx0},${sy0} A ${R},${R} 0 0 ${sweep} ${ex},${ey}`;
        const ux = -Math.sin(aEnd);
        const uy = -Math.cos(aEnd);
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

      let ux = 0;
      let uy = 0;
      if (ld.direction === "x") ux = Math.sign(ld.value);
      else if (ld.direction === "y") uy = -Math.sign(ld.value);
      const arrow = arrowGeom(p.x, p.y, ux, uy, 28);
      return {
        ...base,
        kind: "arrow",
        x: p.x,
        y: p.y,
        ...arrow,
        labelX: p.x + ux * 36 + 6,
        labelY: p.y + uy * 36 - 4,
      };
    })
    .filter(Boolean);
}

/** 生成箭头几何: 起点(sx,sy), 单位方向(fx,fy), 长度 L */
function arrowGeom(sx, sy, fx, fy, L) {
  const nrm = Math.hypot(fx, fy) || 1;
  const dx = fx / nrm;
  const dy = fy / nrm;
  const tx = sx + dx * L;
  const ty = sy + dy * L;
  const px = -dy;
  const py = dx;
  const hx = tx - dx * 11;
  const hy = ty - dy * 11;
  return {
    tx,
    ty,
    arrowPoints: `${tx},${ty} ${hx - px * 5},${hy - py * 5} ${hx + px * 5},${hy + py * 5}`,
  };
}

function fmtLoad(v) {
  if (typeof v !== "number" || Number.isNaN(v)) return String(v ?? 0);
  if (Math.abs(v) >= 100000) return (v / 1000).toFixed(0) + "k";
  if (Math.abs(v) >= 1000) return (v / 1000).toFixed(1) + "k";
  return String(Math.round(v * 100) / 100);
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
