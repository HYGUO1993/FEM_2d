import { useEffect, useRef, useState } from "react";
import { useStore } from "../store.js";
import { screenToWorld, pointSegDist, clamp } from "../canvas/transform.js";
import {
  computeNodes,
  computeElements,
  computeLoads,
  computeConstraints,
  computeDeformed,
  computeForceDiagram,
} from "../canvas/render.js";

const NODE_HIT = 12;
const ELEM_HIT = 8;
const LOAD_HIT = 14;

const TOOL_HINTS = {
  select: "选择: 点击节点/杆件选中, 拖拽节点移动, 空白拖拽平移",
  node: "节点: 点击空白处添加节点 (吸附 0.1m)",
  element: "杆件: 依次点击两个节点创建单元, 点空白取消",
  load: "荷载: 点击节点设置集中力",
  constraint: "约束: 点击节点, 在右侧勾选 ux/uy/rz",
  erase: "删除: 点击荷载/约束/节点/杆件删除",
};

function svgPointOf(svg, e) {
  const rect = svg.getBoundingClientRect();
  return { x: e.clientX - rect.left, y: e.clientY - rect.top };
}

export default function CanvasView() {
  const tool = useStore((s) => s.tool);
  const model = useStore((s) => s.model);
  const view = useStore((s) => s.view);
  const selection = useStore((s) => s.selection);
  const results = useStore((s) => s.results);
  const solved = useStore((s) => s.solved);
  const deformScale = useStore((s) => s.deformScale);
  const pendingNodeA = useStore((s) => s.pendingNodeA);
  const loadDialog = useStore((s) => s.loadDialog);
  // 内力图显示模式: null | "N" | "V" | "M"
  const [diagramMode, setDiagramMode] = useState(null);

  const wrapRef = useRef(null);
  const svgRef = useRef(null);
  const drag = useRef({ mode: "none" });

  // 首屏自适应 + 尺寸观察
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const st = useStore.getState();
    st.setViewBox(el.clientWidth, el.clientHeight);
    st.fitView(el.clientWidth, el.clientHeight);
    const ro = new ResizeObserver(() => {
      useStore.getState().setViewBox(el.clientWidth, el.clientHeight);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  // 滚轮缩放（被动监听, 以鼠标为锚点）
  useEffect(() => {
    const el = svgRef.current;
    if (!el) return;
    const handler = (e) => {
      e.preventDefault();
      const st = useStore.getState();
      const v = st.view;
      if (!v) return;
      const { x: mx, y: my } = svgPointOf(el, e);
      const factor = e.deltaY < 0 ? 1.15 : 1 / 1.15;
      const ns = clamp(v.scale * factor, 2, 2000);
      const k = ns / v.scale;
      st.setView({
        scale: ns,
        panX: mx - (mx - v.panX) * k,
        panY: my - (my - v.panY) * k,
      });
    };
    el.addEventListener("wheel", handler, { passive: false });
    return () => el.removeEventListener("wheel", handler);
  }, []);

  // —— 命中检测 ——
  function hitNode(sx, sy) {
    const nodes = computeNodes(model, view, selection);
    for (let i = nodes.length - 1; i >= 0; i--) {
      const n = nodes[i];
      if (Math.hypot(sx - n.cx, sy - n.cy) <= NODE_HIT) return n;
    }
    return null;
  }
  function hitElement(sx, sy) {
    const els = computeElements(model, view, selection);
    for (const e of els) {
      if (pointSegDist(sx, sy, e.x1, e.y1, e.x2, e.y2) <= ELEM_HIT) return e;
    }
    return null;
  }
  function hitLoad(sx, sy) {
    const loads = computeLoads(model, view, selection);
    for (const l of loads) {
      const tx = l.kind === "rz" ? l.head.x : l.tx;
      const ty = l.kind === "rz" ? l.head.y : l.ty;
      if (Math.hypot(sx - tx, sy - ty) <= LOAD_HIT) return l;
      if (Math.hypot(sx - l.x, sy - l.y) <= LOAD_HIT) return l;
    }
    return null;
  }
  function hitConstraint(sx, sy) {
    const cons = computeConstraints(model, view, selection);
    for (const c of cons) {
      if (Math.abs(sx - c.x) <= 12 && sy >= c.y - 6 && sy <= c.y + c.h + 10) return c;
    }
    return null;
  }

  function onMouseDown(e) {
    const st = useStore.getState();
    const v = st.view;
    if (!v) return;
    const { x: sx, y: sy } = svgPointOf(svgRef.current, e);

    // 中键任意模式下平移
    if (e.button === 1) {
      drag.current = { mode: "pan", sx, sy, startView: { ...v } };
      return;
    }
    if (e.button !== 0) return;

    const curTool = st.tool;

    if (curTool === "erase") {
      const ld = hitLoad(sx, sy);
      if (ld) return st.deleteLoad(ld.id);
      const c = hitConstraint(sx, sy);
      if (c) return st.removeConstraint(c.node);
      const n = hitNode(sx, sy);
      if (n) return st.deleteNode(n.id);
      const el = hitElement(sx, sy);
      if (el) return st.deleteElement(el.id);
      return;
    }

    if (curTool === "node") {
      if (hitNode(sx, sy)) return; // 避免叠点
      const w = screenToWorld(sx, sy, v);
      const x = Math.round(w.x / 0.1) * 0.1;
      const y = Math.round(w.y / 0.1) * 0.1;
      st.addNode(Number(x.toFixed(4)), Number(y.toFixed(4)));
      return;
    }

    if (curTool === "element") {
      const n = hitNode(sx, sy);
      if (n) {
        if (st.pendingNodeA == null) st.setPendingNodeA(n.id);
        else if (st.pendingNodeA !== n.id) {
          st.addElement(st.pendingNodeA, n.id);
          st.setPendingNodeA(null);
        }
      } else {
        st.setPendingNodeA(null); // 点空白取消
      }
      return;
    }

    if (curTool === "load") {
      const n = hitNode(sx, sy);
      if (n) st.setLoadDialog({ node: n.id });
      return;
    }

    if (curTool === "constraint") {
      const n = hitNode(sx, sy);
      if (n) st.select({ type: "constraint", id: n.id });
      return;
    }

    // select
    const n = hitNode(sx, sy);
    if (n) {
      st.select({ type: "node", id: n.id });
      drag.current = { mode: "node", nodeId: n.id, sx, sy };
      return;
    }
    const el = hitElement(sx, sy);
    if (el) {
      st.select({ type: "element", id: el.id });
      return;
    }
    st.select(null);
    drag.current = { mode: "pan", sx, sy, startView: { ...v } };
  }

  function onMouseMove(e) {
    const d = drag.current;
    const st = useStore.getState();
    const v = st.view;
    if (!v || d.mode === "none") return;
    const { x: sx, y: sy } = svgPointOf(svgRef.current, e);
    if (d.mode === "pan") {
      st.setView({
        ...d.startView,
        panX: d.startView.panX + (sx - d.sx),
        panY: d.startView.panY + (sy - d.sy),
      });
    } else if (d.mode === "node") {
      const w = screenToWorld(sx, sy, v);
      const x = Math.round(w.x / 0.05) * 0.05;
      const y = Math.round(w.y / 0.05) * 0.05;
      st.moveNode(d.nodeId, Number(x.toFixed(4)), Number(y.toFixed(4)));
    }
  }

  function endDrag() {
    drag.current = { mode: "none" };
  }

  const nodes = computeNodes(model, view, selection);
  const elements = computeElements(model, view, selection);
  const loads = computeLoads(model, view, selection);
  const constraints = computeConstraints(model, view, selection);
  const deformed =
    solved && results ? computeDeformed(model, results, deformScale, view) : [];
  const forceDiagrams =
    solved && results && diagramMode ? computeForceDiagram(model, results, diagramMode, view) : [];

  const pendingNode = pendingNodeA != null ? nodes.find((n) => n.id === pendingNodeA) : null;

  const DIAGRAM_META = {
    N: { label: "轴力图 N", color: "#4da3ff", unit: "N" },
    V: { label: "剪力图 V", color: "#fb923c", unit: "N" },
    M: { label: "弯矩图 M", color: "#f87171", unit: "N·m" },
  };

  return (
    <div className="canvas-wrap" ref={wrapRef}>
      <svg
        ref={svgRef}
        className="canvas-svg"
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={endDrag}
        onMouseLeave={endDrag}
        onContextMenu={(e) => e.preventDefault()}
      >
        {/* 变形图（红色虚线, 最底层） */}
        {deformed.map((d, i) => (
          <line key={`def-${i}`} className="deformed" x1={d.x1} y1={d.y1} x2={d.x2} y2={d.y2} />
        ))}

        {/* 内力图（N/V/M, 叠加在杆件下层） */}
        {forceDiagrams.map((d) => (
          <path
            key={d.key}
            className={`force-diagram ${diagramMode}`}
            d={d.d}
            fill="none"
            stroke={DIAGRAM_META[diagramMode]?.color}
          />
        ))}

        {/* 杆件 */}
        {elements.map((e) => (
          <line
            key={`el-${e.id}`}
            className={`element ${e.type === "truss" ? "truss" : "frame"} ${
              e.selected ? "selected" : ""
            }`}
            x1={e.x1}
            y1={e.y1}
            x2={e.x2}
            y2={e.y2}
          />
        ))}

        {/* 荷载 */}
        {loads.map((l) => (
          <g key={`ld-${l.id}`} className={l.selected ? "selected" : ""}>
            {l.kind === "rz" ? (
              <path className="load-line" d={l.d} fill="none" />
            ) : (
              <line className="load-line" x1={l.x} y1={l.y} x2={l.tx} y2={l.ty} />
            )}
            <polygon className="load-head" points={l.arrowPoints} />
            <text className="load-label" x={l.labelX ?? l.head.x + 8} y={l.labelY ?? l.head.y - 6}>
              {l.value}
            </text>
          </g>
        ))}

        {/* 约束 */}
        {constraints.map((c) => (
          <g key={`c-${c.node}`} className={c.selected ? "selected" : ""}>
            {c.hasRz && (
              <circle className="constraint-rz" cx={c.x} cy={c.y} r={7} fill="none" />
            )}
            <polygon
              className={`constraint ${c.fixed ? "fixed" : "roller"}`}
              points={c.pts}
            />
            {c.roller && (
              <circle className="constraint-roller-dot" cx={c.x} cy={c.y + c.h} r={3} />
            )}
          </g>
        ))}

        {/* 节点 */}
        {nodes.map((n) => (
          <g key={`n-${n.id}`}>
            <circle
              className={`node ${n.type === "truss" ? "truss" : "frame"} ${
                n.selected ? "selected" : ""
              }`}
              cx={n.cx}
              cy={n.cy}
              r={5}
            />
            <text className="node-label" x={n.cx + 7} y={n.cy - 7}>
              {n.id}
            </text>
          </g>
        ))}

        {/* element 模式第一个点高亮 */}
        {pendingNode && (
          <circle className="pending-node" cx={pendingNode.cx} cy={pendingNode.cy} r={10} fill="none" />
        )}
      </svg>

      {loadDialog && <LoadDialog />}

      <div className="canvas-hint">{TOOL_HINTS[tool]}</div>

      {solved && results && (
        <div className="diagram-bar">
          <button
            className={`btn small ${diagramMode === null ? "active" : ""}`}
            onClick={() => setDiagramMode(null)}
            title="隐藏内力图"
          >
            原图
          </button>
          <button
            className={`btn small ${diagramMode === "N" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "N" ? null : "N")}
            title="轴力图 N"
          >
            轴力 N
          </button>
          <button
            className={`btn small ${diagramMode === "V" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "V" ? null : "V")}
            title="剪力图 V"
          >
            剪力 V
          </button>
          <button
            className={`btn small ${diagramMode === "M" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "M" ? null : "M")}
            title="弯矩图 M"
          >
            弯矩 M
          </button>
        </div>
      )}
    </div>
  );
}

// —— 荷载内联对话框 ——
function LoadDialog() {
  const [direction, setDirection] = useState("y");
  const [value, setValue] = useState("-10000");
  const dialog = useStore((s) => s.loadDialog);

  if (!dialog) return null;
  const st = useStore.getState();

  return (
    <div className="load-dialog">
      <div className="load-dialog-title">集中力 @节点 {dialog.node}</div>
      <label className="field">
        <span>方向</span>
        <select value={direction} onChange={(e) => setDirection(e.target.value)}>
          <option value="x">x (水平)</option>
          <option value="y">y (竖向)</option>
          <option value="rz">rz (弯矩)</option>
        </select>
      </label>
      <label className="field">
        <span>数值</span>
        <input
          type="number"
          value={value}
          onChange={(e) => setValue(e.target.value)}
        />
      </label>
      <div className="load-dialog-actions">
        <button
          className="btn primary small"
          onClick={() => {
            const v = Number(value);
            if (Number.isNaN(v)) return;
            st.addLoad(dialog.node, direction, v);
          }}
        >
          确认
        </button>
        <button className="btn small" onClick={() => st.setLoadDialog(null)}>
          取消
        </button>
      </div>
    </div>
  );
}
