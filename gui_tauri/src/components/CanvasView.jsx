import { useEffect, useRef, useState } from "react";
import { useStore } from "../store.js";
import { screenToWorld, pointSegDist, clamp } from "../canvas/transform.js";
import { LOAD_TYPES } from "../model.js";
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
      else {
        const el = hitElement(sx, sy);
        if (el) st.setLoadDialog({ node: -1, element: el.id });
      }
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

  // 节点位移标注: 求解后显示 ux/uy (mm 或 m), 若变形图可见
  const displacementMap = new Map(
    (results?.displacements || []).map((d) => [d.node, d])
  );

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
            {l.kind === "rz" && (
              <>
                <path className="load-line" d={l.d} fill="none" />
                <polygon className="load-head" points={l.arrowPoints} />
              </>
            )}
            {l.kind === "arrow" && (
              <>
                <line className="load-line" x1={l.x} y1={l.y} x2={l.tx} y2={l.ty} />
                <polygon className="load-head" points={l.arrowPoints} />
              </>
            )}
            {l.kind === "elemArrow" && (
              <>
                <line className="load-line" x1={l.sx} y1={l.sy} x2={l.tx} y2={l.ty} />
                <polygon className="load-head" points={l.arrowPoints} />
              </>
            )}
            {l.kind === "udl" && (
              <>
                {l.arrows.map((a, i) => (
                  <g key={i}>
                    <line className="load-line" x1={a.sx} y1={a.sy} x2={a.tx} y2={a.ty} />
                    <polygon className="load-head" points={a.arrowPoints} />
                  </g>
                ))}
              </>
            )}
            {l.kind === "temp" && (
              <text className="temp-label" x={l.x} y={l.y}>
                {l.label}
              </text>
            )}
            {l.labelX != null && (
              <text className="load-label" x={l.labelX} y={l.labelY}>
                {l.label ?? l.value}
              </text>
            )}
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

        {/* 节点位移标注 (求解后, 画在节点下方) */}
        {solved && results && !diagramMode && (
          <g className="disp-labels">
            {nodes.map((n) => {
              const d = displacementMap.get(n.id);
              if (!d) return null;
              const ux = d.ux || 0;
              const uy = d.uy || 0;
              if (Math.abs(ux) < 1e-12 && Math.abs(uy) < 1e-12) return null;
              const mag = Math.hypot(ux, uy);
              const mm = mag >= 1e-3; // >1mm 用 mm 显示
              const txt = mm
                ? `${(mag * 1000).toFixed(1)}mm`
                : `${(mag * 1000).toFixed(3)}mm`;
              return (
                <text key={`dl-${n.id}`} className="disp-label" x={n.cx} y={n.cy + 16}>
                  {txt}
                </text>
              );
            })}
          </g>
        )}

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

// —— 荷载对话框(支持完整荷载类型) ——
function LoadDialog() {
  const [type, setType] = useState("nodalForce");
  const [direction, setDirection] = useState("y");
  const [value, setValue] = useState("-10000");
  const [element, setElement] = useState("0");
  const [position, setPosition] = useState("1");
  const [T0, setT0] = useState("0");
  const [T1, setT1] = useState("0");
  const dialog = useStore((s) => s.loadDialog);
  const model = useStore((s) => s.model);

  if (!dialog) return null;
  const st = useStore.getState();

  const isNodeLoad = type === "nodalForce" || type === "momentOnPoint";
  const isElementLoad = !isNodeLoad; // 单元荷载(含温度)都需要选单元
  const isTemp = type === "temperature";

  function submit() {
    const v = Number(value);
    if (Number.isNaN(v)) return;
    const params = { type, direction, value: v };
    if (isNodeLoad) {
      params.node = dialog.node;
    } else if (isElementLoad) {
      const eid = Number(element);
      if (Number.isNaN(eid) || !model.elements.some((e) => e.id === eid)) {
        st.setToast("请选择有效的单元", true);
        return;
      }
      params.element = eid;
      const pos = Number(position);
      if (!Number.isNaN(pos)) params.position = pos;
      // 节点荷载(横向集中力等)需要节点 = 单元起点
      if (type === "lateralForce" || type === "axialForce") {
        const el = model.elements.find((e) => e.id === eid);
        params.node = el ? el.nodeI : -1;
      }
    } else if (isTemp) {
      params.T0 = Number(T0) || 0;
      params.T1 = Number(T1) || 0;
      params.node = -1;
      // 温度荷载必须指定单元 (后端按单元算固端力)
      const eid = Number(element);
      if (Number.isNaN(eid) || !model.elements.some((e) => e.id === eid)) {
        st.setToast("温度荷载需要选择作用单元", true);
        return;
      }
      params.element = eid;
    }
    st.addLoad(params);
  }

  return (
    <div className="load-dialog">
      <div className="load-dialog-title">
        添加荷载 {isNodeLoad ? `@节点 ${dialog.node}` : ""}
      </div>
      <label className="field">
        <span>类型</span>
        <select value={type} onChange={(e) => setType(e.target.value)}>
          {Object.entries(LOAD_TYPES).map(([k, label]) => (
            <option key={k} value={k}>
              {label}
            </option>
          ))}
        </select>
      </label>
      {!isTemp && (
        <label className="field">
          <span>方向</span>
          <select value={direction} onChange={(e) => setDirection(e.target.value)}>
            <option value="x">x (水平)</option>
            <option value="y">y (竖向)</option>
            <option value="rz">rz (转动)</option>
          </select>
        </label>
      )}
      {isElementLoad && (
        <label className="field">
          <span>作用单元</span>
          <select value={element} onChange={(e) => setElement(e.target.value)}>
            {model.elements.map((e) => (
              <option key={e.id} value={e.id}>
                单元 #{e.id} (节点 {e.nodeI}–{e.nodeJ})
              </option>
            ))}
          </select>
        </label>
      )}
      {isTemp ? (
        <>
          <label className="field">
            <span>下表面温变 T0 (℃)</span>
            <input type="number" value={T0} onChange={(e) => setT0(e.target.value)} />
          </label>
          <label className="field">
            <span>上表面温变 T1 (℃)</span>
            <input type="number" value={T1} onChange={(e) => setT1(e.target.value)} />
          </label>
        </>
      ) : (
        <>
          <label className="field">
            <span>数值 {isElementLoad ? "(N/m)" : "(N 或 N·m)"}</span>
            <input type="number" value={value} onChange={(e) => setValue(e.target.value)} />
          </label>
          {isElementLoad && (
            <label className="field">
              <span>作用长度/位置 position (m, 从单元起点, 默认整跨)</span>
              <input type="number" value={position} onChange={(e) => setPosition(e.target.value)} />
            </label>
          )}
        </>
      )}
      <div className="load-dialog-actions">
        <button className="btn primary small" onClick={submit}>
          确认
        </button>
        <button className="btn small" onClick={() => st.setLoadDialog(null)}>
          取消
        </button>
      </div>
    </div>
  );
}
