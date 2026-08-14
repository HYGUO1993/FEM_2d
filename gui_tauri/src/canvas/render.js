// 模型 → SVG 元素几何计算（GUI_REDESIGN_PLAN §4.6）
// 返回描述对象数组, 由 CanvasView 渲染成 JSX/SVG 元素。
// 注意: 本模块不依赖 i18n —— label 字段返回语义键(labelKey), 由 CanvasView 用 t() 翻译。

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
          labelKey: "canvas.tempLabel",
          labelArgs: { a: ld.T0 ?? 0, b: ld.T1 ?? 0 },
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
            labelKey: ld.type === "axialPressure" ? "canvas.axialLabel" : "canvas.udlLabel",
            labelArgs: { v: fmtLoad(ld.value) },
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
            labelKey: "canvas.linearLabel",
            labelArgs: { v: fmtLoad(ld.value) },
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

      // 支座位移: 蓝绿色位移箭头
      if (ld.type === "supportMove") {
        let ux = 0;
        let uy = 0;
        if (ld.direction === "x") ux = Math.sign(ld.value);
        else if (ld.direction === "y") uy = -Math.sign(ld.value);
        else if (ld.direction === "rz") { /* 转动位移: 圆弧 */ }
        const g = arrowGeom(p.x, p.y - 20, ux, uy, 20);
        return {
          ...base,
          kind: "supportMove",
          x: p.x,
          y: p.y - 20,
          ...g,
          label: `${ld.value} m`,
          labelX: p.x + ux * 26 + 6,
          labelY: p.y - 20 + uy * 26 - 4,
        };
      }

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
      const hasUx = c.dofs.includes("ux");
      const hasUy = c.dofs.includes("uy");
      const hasRz = c.dofs.includes("rz");

      // 结构力学惯例分类
      let kind;
      if (hasUx && hasUy && hasRz) kind = "fixed";       // 固定端(固支)
      else if (hasUx && hasUy) kind = "pinned";           // 固定铰支座
      else if (hasUy && !hasUx) kind = "roller";          // 滚轴支座(活动铰)
      else if (hasUx && !hasUy) kind = "rollerX";         // 水平滚轴
      else if (hasRz && !hasUx && !hasUy) kind = "rzOnly"; // 仅转角约束
      else kind = "partial";                              // 其它组合

      const h = 11;      // 支座高度
      const half = 9;    // 半宽
      const y = p.y + 9; // 支座顶部(紧贴节点下方)

      // —— 固定端(固支): 垂直杆件一条横线, 横线下 45° 同方向平行斜线 ——
      const barY = y;                    // 横线(与杆件连接处)
      const barX1 = p.x - half;
      const barX2 = p.x + half;
      // 斜线: 45° 向右下, 均匀分布, 同方向平行
      const slantH = 8;                  // 斜线垂直高度(=水平跨度, 45°)
      const slantN = 4;                  // 斜线数量
      const slantLines = [];
      for (let i = 0; i < slantN; i++) {
        const sx = barX1 + (i * (barX2 - barX1 - slantH)) / (slantN - 1);
        slantLines.push([sx, barY, sx + slantH, barY + slantH]);
      }

      // —— 铰链符号: 两边圆圈中间细线相连 ——
      const pinY = y + 6;                // 铰链中心
      const pinR = 5;                    // 圆半径
      const pinGap = 8;                  // 两圆间距
      const pinX1 = p.x - pinGap;
      const pinX2 = p.x + pinGap;
      const linkY = y + 2;               // 上部连接线(与杆件)

      return {
        node: c.node,
        dofs: c.dofs,
        x: p.x,
        y,
        kind,
        hasUx,
        hasUy,
        hasRz,
        // 固定端几何
        barY, barX1, barX2,
        slantH, slantLines,
        // 铰链几何
        pinY, pinR, pinX1, pinX2, linkY,
        // 地面线 (铰支/滚轴下方)
        groundY: pinY + pinR + 4,
        // 滚轮 (滚轴支座): 两个小圆
        rollerY: pinY + pinR + 4,
        // rz 转角约束圆弧
        rzArc: hasRz ? `M ${p.x - 8},${p.y + 2} A 8,8 0 0 1 ${p.x + 8},${p.y + 2}` : null,
        selected: selection && selection.type === "constraint" && selection.id === c.node,
      };
    })
    .filter(Boolean);
}

/**
 * 变形图: 沿杆件用 Hermite 三次插值 (利用节点位移 ux/uy + 转角 rz),
 * 画出真实弯曲形状 (直线连接会丢失均布荷载下的挠曲形态)。
 * 局部坐标: x' 沿杆轴, y' 左侧法向; rz = dv/dx' (小变形)。
 * 返回 {key, d} 路径, 由 CanvasView 渲染为 <path>。
 */
export function computeDeformed(model, results, deformScale, view) {
  if (!results || !Array.isArray(results.displacements)) return [];
  const disp = new Map(results.displacements.map((d) => [d.node, d]));
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  const SEG = 12;
  const out = [];
  for (const e of model.elements || []) {
    const a = byId.get(e.nodeI);
    const b = byId.get(e.nodeJ);
    if (!a || !b) continue;
    const da = disp.get(e.nodeI);
    const db = disp.get(e.nodeJ);
    const uxa = da ? (da.ux || 0) * deformScale : 0;
    const uya = da ? (da.uy || 0) * deformScale : 0;
    const rza = da ? da.rz || 0 : 0;
    const uxb = db ? (db.ux || 0) * deformScale : 0;
    const uyb = db ? (db.uy || 0) * deformScale : 0;
    const rzb = db ? db.rz || 0 : 0;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.hypot(dx, dy) || 1e-9;
    const c = dx / len;
    const s = dy / len;
    // 局部位移 (i/j 端): 轴向 u、横向 v、斜率 dv/dt = rz*deformScale*len
    const ua = c * uxa + s * uya;
    const va = -s * uxa + c * uya;
    const ub = c * uxb + s * uyb;
    const vb = -s * uxb + c * uyb;
    const ta = rza * deformScale * len;
    const tb = rzb * deformScale * len;
    const pts = [];
    for (let k = 0; k <= SEG; k++) {
      const t = k / SEG;
      const t2 = t * t;
      const t3 = t2 * t;
      const h00 = 2 * t3 - 3 * t2 + 1;
      const h10 = t3 - 2 * t2 + t;
      const h01 = -2 * t3 + 3 * t2;
      const h11 = t3 - t2;
      const ux = ua + (ub - ua) * t; // 轴向线性
      const v = h00 * va + h10 * ta + h01 * vb + h11 * tb; // 横向 Hermite
      const wx = a.x + (t * len + ux) * c + v * -s;
      const wy = a.y + (t * len + ux) * s + v * c;
      const p = worldToScreen(wx, wy, view);
      pts.push(`${p.x},${p.y}`);
    }
    out.push({ key: `def-${e.id}`, d: `M ${pts.join(" L ")}` });
  }
  return out;
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
 * 内力图: 按单元荷载精确绘制杆内 N/V/M 分布 (不再线性直连两端)。
 * kind: "N" | "V" | "M"
 *
 * 符号约定（与求解器端力一致, 局部坐标 x' 沿杆轴, y' 左侧法向）:
 *  - M: i 端取 M_i, 杆内 M(x) = M_i - V_i·x - Σ(荷载弯矩贡献), j 端自动等于 -M_j
 *    (代数恒等: M(L) = -M_j), 跨节点连续
 *  - V: V(x) = V_i + Σ(荷载剪力贡献), j 端自动等于 -V_j
 *  - N: N(x) = N_i + Σ(轴向荷载贡献)
 * 荷载贡献 (value 为局部坐标带符号值, 与 FixedEndForceCalcu 同约定):
 *  - lateralUniformPressure q: V += q·x, M += q·x²/2  (x 截断于作用长度 a)
 *  - lateralForce P@a:        V += P·[x≥a], M += P·(x-a)·[x≥a]
 *  - lateralLinearlyPressure 0→q: V += q·x²/2a, M += q·x³/6a
 *  - axialPressure q:         N += q·x
 *  - axialForce P@a:          N += P·[x≥a]
 * 正弯矩画在杆轴上方(受拉侧), 负画下方; V/N 正负分侧。
 */
export function loadContribAt(loads, len, x) {
  let V = 0;
  let M = 0;
  let N = 0;
  for (const ld of loads) {
    const a = ld.position || 0;
    if (ld.type === "lateralUniformPressure") {
      const aa = a <= 0 || a >= len ? len : Math.min(a, len);
      const xx = Math.min(x, aa);
      V += ld.value * xx;
      M += (ld.value * xx * xx) / 2;
    } else if (ld.type === "lateralForce") {
      if (x >= a) {
        V += ld.value;
        M += ld.value * (x - a);
      }
    } else if (ld.type === "lateralLinearlyPressure") {
      const aa = a <= 0 || a >= len ? len : Math.min(a, len);
      const xx = Math.min(x, aa);
      V += (ld.value * xx * xx) / (2 * aa);
      M += (ld.value * xx * xx * xx) / (6 * aa);
    } else if (ld.type === "axialPressure") {
      const aa = a <= 0 || a >= len ? len : Math.min(a, len);
      N += ld.value * Math.min(x, aa);
    } else if (ld.type === "axialForce") {
      if (x >= a) N += ld.value;
    }
  }
  return { V, M, N };
}

export function computeForceDiagram(model, results, kind, view) {
  if (!results || !Array.isArray(results.endForces)) return [];
  const byId = new Map((model.nodes || []).map((n) => [n.id, n]));
  const forceById = new Map(results.endForces.map((f) => [f.element, f]));
  const loadsByElem = new Map();
  for (const ld of model.loads || []) {
    if (ld.element == null || ld.element < 0) continue;
    if (!loadsByElem.has(ld.element)) loadsByElem.set(ld.element, []);
    loadsByElem.get(ld.element).push(ld);
  }

  // —— 1) 逐单元采样内力分布 ——
  const segs = [];
  let maxVal = 1e-9;
  for (const el of model.elements || []) {
    const a = byId.get(el.nodeI);
    const b = byId.get(el.nodeJ);
    const f = forceById.get(el.id);
    if (!a || !b || !f) continue;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.hypot(dx, dy) || 1e-9;
    const loads = loadsByElem.get(el.id) || [];
    const ni = f.nodeI || {};
    const nj = f.nodeJ || {};
    const Mi = ni.M || 0;
    const Vi = ni.V || 0;
    const Ni = ni.N || 0;

    // 采样位置: 端点 + 集中荷载位置(剪/轴力跳变) + 均匀细分
    const xs = [0, len];
    for (const ld of loads) {
      if ((ld.type === "lateralForce" || ld.type === "axialForce") && ld.position > 0 && ld.position < len) {
        xs.push(ld.position);
      }
    }
    const SEG = 16;
    for (let k = 1; k < SEG; k++) xs.push((len * k) / SEG);
    xs.sort((p, q) => p - q);

    const pts = [];
    let lastX = -1;
    for (const x of xs) {
      if (x === lastX) continue;
      lastX = x;
      const cc = loadContribAt(loads, len, x);
      let v;
      if (kind === "M") v = Mi - Vi * x - cc.M;
      else if (kind === "V") v = Vi + cc.V;
      else v = Ni + cc.N;
      maxVal = Math.max(maxVal, Math.abs(v));
      pts.push({ x, v });
    }
    segs.push({ el, a, b, len, pts });
  }

  // —— 2) 归一化并映射到屏幕 ——
  const out = [];
  for (const seg of segs) {
    const { el, a, b, len, pts } = seg;
    const ux = (b.x - a.x) / len;
    const uy = (b.y - a.y) / len;
    const nx = -uy; // 杆轴法向 (世界), 正偏移 → 杆轴一侧
    const ny = ux;
    const k = (0.25 * len) / maxVal;
    const d =
      "M " +
      pts
        .map((p) => {
          const s = worldToScreen(a.x + nx * p.v * k, a.y + ny * p.v * k, view);
          return `${s.x},${s.y}`;
        })
        .join(" L ");
    out.push({ key: `fg-${kind}-${el.id}`, d, maxVal, samples: pts });
  }
  return out;
}
