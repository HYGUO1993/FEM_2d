import { useEffect, useMemo, useRef, useState } from "react";
import { useStore } from "../store.js";
import { screenToWorld, pointSegDist, clamp } from "../canvas/transform.js";
import { LOAD_TYPES } from "../model.js";
import { t } from "../i18n.js";
import * as ipc from "../ipc.js";
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
  select: "hint.select",
  node: "hint.node",
  element: "hint.element",
  load: "hint.load",
  constraint: "hint.constraint",
  erase: "hint.erase",
};

const DIAGRAM_META = {
  N: { color: "#4da3ff" },
  V: { color: "#fb923c" },
  M: { color: "#f87171" },
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
  // 显示选项
  const displayOptions = useStore((s) => s.displayOptions);

  const wrapRef = useRef(null);
  const svgRef = useRef(null);
  const drag = useRef({ mode: "none" });

  // 导出画布为 PNG (SVG → canvas → 下载)
  function exportPng() {
    const svg = svgRef.current;
    if (!svg) return;
    try {
      const rect = svg.getBoundingClientRect();
      const w = rect.width;
      const h = rect.height;
      // 克隆 SVG, 内联样式与背景
      const clone = svg.cloneNode(true);
      clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      clone.setAttribute("width", w);
      clone.setAttribute("height", h);
      // 计算样式背景色 (读 CSS 变量)
      const bg = getComputedStyle(document.documentElement).getPropertyValue("--bg").trim() || "#0f1115";
      const svgStr =
        `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">` +
        `<rect width="${w}" height="${h}" fill="${bg}"/>` +
        clone.outerHTML +
        `</svg>`;
      const blob = new Blob([svgStr], { type: "image/svg+xml;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const img = new Image();
      img.onload = () => {
        const canvas = document.createElement("canvas");
        canvas.width = w * 2; // 2x 高清
        canvas.height = h * 2;
        const ctx = canvas.getContext("2d");
        ctx.scale(2, 2);
        ctx.drawImage(img, 0, 0, w, h);
        URL.revokeObjectURL(url);
        // 弹原生保存对话框, 由后端写入文件 (a.click 在 WebView2 不可靠)
        const b64 = canvas.toDataURL("image/png");
        ipc.exportPngFile(`femlab_${Date.now()}`, b64)
          .then((path) => useStore.getState().setToast(t("canvas.toast.exported") + " → " + path))
          .catch((e) => {
            const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
            if (msg !== "用户取消了保存") {
              useStore.getState().setToast(t("canvas.toast.exportErr", { msg }), true);
            }
          });
      };
      img.onerror = () => useStore.getState().setToast(t("canvas.toast.exportFailed"), true);
      img.src = url;
    } catch (e) {
      useStore.getState().setToast(t("canvas.toast.exportErr", { msg: e.message }), true);
    }
  }

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
      if (n) st.setLoadDialog({ node: n.id, element: -1 });
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

  const nodes = useMemo(() => computeNodes(model, view, selection), [model, view, selection]);
  const elements = useMemo(() => computeElements(model, view, selection), [model, view, selection]);
  const loads = useMemo(() => computeLoads(model, view, selection), [model, view, selection]);
  const constraints = useMemo(() => computeConstraints(model, view, selection), [model, view, selection]);
  const deformed = useMemo(
    () => (solved && results ? computeDeformed(model, results, deformScale, view) : []),
    [solved, results, model, deformScale, view]
  );
  const forceDiagrams = useMemo(
    () => (solved && results && diagramMode ? computeForceDiagram(model, results, diagramMode, view) : []),
    [solved, results, model, diagramMode, view]
  );

  // 节点位移标注: 求解后显示 ux/uy (mm 或 m), 若变形图可见
  const displacementMap = useMemo(
    () => new Map((results?.displacements || []).map((d) => [d.node, d])),
    [results]
  );

  const pendingNode = pendingNodeA != null ? nodes.find((n) => n.id === pendingNodeA) : null;

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
        {/* 变形图（红色虚线, 最底层, Hermite 三次插值曲线） */}
        {deformed.map((d) => (
          <path key={d.key} className="deformed" d={d.d} fill="none" />
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
          <g key={`el-${e.id}`}>
            <line
              className={`element ${e.type === "truss" ? "truss" : "frame"} ${
                e.selected ? "selected" : ""
              }`}
              x1={e.x1}
              y1={e.y1}
              x2={e.x2}
              y2={e.y2}
            />
            {displayOptions.elementLabels && (
              <text
                className="elem-label"
                x={(e.x1 + e.x2) / 2}
                y={(e.y1 + e.y2) / 2 - 6}
              >
                {e.id}
              </text>
            )}
          </g>
        ))}

        {/* 荷载 */}
        {displayOptions.loads &&
          loads.map((l) => (
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
            {l.kind === "supportMove" && (
              <>
                <line className="support-move-line" x1={l.x} y1={l.y} x2={l.tx} y2={l.ty} />
                <polygon className="support-move-head" points={l.arrowPoints} />
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
                {l.labelKey ? t(l.labelKey, l.labelArgs) : l.label}
              </text>
            )}
            {l.labelX != null && (
              <text className="load-label" x={l.labelX} y={l.labelY}>
                {l.labelKey ? t(l.labelKey, l.labelArgs) : (l.label ?? l.value)}
              </text>
            )}
          </g>
        ))}

        {/* 约束 */}
        {displayOptions.constraints &&
          constraints.map((c) => (
          <g key={`c-${c.node}`} className={c.selected ? "selected" : ""}>
            {/* 转角约束圆弧 */}
            {c.rzArc && <path className="constraint-rz" d={c.rzArc} fill="none" />}

            {c.kind === "fixed" && (
              <>
                {/* 固定端: 横线 + 45° 同方向平行斜线 (教材标准) */}
                <line className="constraint-bar" x1={c.barX1} y1={c.barY} x2={c.barX2} y2={c.barY} />
                {c.slantLines.map(([x1, y1, x2, y2], i) => (
                  <line key={i} className="constraint-slant" x1={x1} y1={y1} x2={x2} y2={y2} />
                ))}
              </>
            )}

            {(c.kind === "pinned" || c.kind === "roller" || c.kind === "rollerX" || c.kind === "partial") && (
              <>
                {/* 铰链: 上部连接线 + 两圆 + 中间细线 */}
                <line className="constraint-bar" x1={c.pinX1} y1={c.linkY} x2={c.pinX2} y2={c.linkY} />
                <circle className="constraint-pin" cx={c.pinX1} cy={c.pinY} r={c.pinR} fill="none" />
                <circle className="constraint-pin" cx={c.pinX2} cy={c.pinY} r={c.pinR} fill="none" />
                <line className="constraint-bar" x1={c.pinX1} y1={c.pinY} x2={c.pinX2} y2={c.pinY} />
                {/* 地面线 */}
                <line className="constraint-ground" x1={c.pinX1 - c.pinR - 4} y1={c.groundY} x2={c.pinX2 + c.pinR + 4} y2={c.groundY} />
                {/* 滚轮 (活动铰): 两圆下方再画滚轮圆 */}
                {c.kind === "roller" && (
                  <>
                    <circle className="constraint-roller-dot" cx={c.pinX1} cy={c.rollerY} r={3} />
                    <circle className="constraint-roller-dot" cx={c.pinX2} cy={c.rollerY} r={3} />
                  </>
                )}
              </>
            )}

            {c.kind === "rzOnly" && (
              <>
                <line className="constraint-ground" x1={c.x - 12} y1={c.groundY} x2={c.x + 12} y2={c.groundY} />
                <circle className="constraint-roller-dot" cx={c.x} cy={c.groundY} r={3} />
              </>
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
            {displayOptions.nodeLabels && (
              <text className="node-label" x={n.cx + 7} y={n.cy - 7}>
                {n.id}
              </text>
            )}
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

      {loadDialog && <LoadDialog key={`${loadDialog.node}-${loadDialog.element}`} />}

      {/* 导出图片按钮 */}
      <button
        className="btn small export-btn"
        onClick={exportPng}
        title={t("canvas.exportPngTip")}
      >
        {t("canvas.exportPng")}
      </button>

      <div className="canvas-hint">{t(TOOL_HINTS[tool])}</div>

      {solved && results && (
        <div className="diagram-bar">
          <button
            className={`btn small ${diagramMode === null ? "active" : ""}`}
            onClick={() => setDiagramMode(null)}
            title={t("canvas.originalTip")}
          >
            {t("canvas.original")}
          </button>
          <button
            className={`btn small ${diagramMode === "N" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "N" ? null : "N")}
            title={t("canvas.nDiagramTip")}
          >
            {t("canvas.nDiagram")}
          </button>
          <button
            className={`btn small ${diagramMode === "V" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "V" ? null : "V")}
            title={t("canvas.vDiagramTip")}
          >
            {t("canvas.vDiagram")}
          </button>
          <button
            className={`btn small ${diagramMode === "M" ? "active" : ""}`}
            onClick={() => setDiagramMode(diagramMode === "M" ? null : "M")}
            title={t("canvas.mDiagramTip")}
          >
            {t("canvas.mDiagram")}
          </button>
        </div>
      )}
    </div>
  );
}

// —— 荷载对话框(支持完整荷载类型) ——
function LoadDialog() {
  const dialog = useStore((s) => s.loadDialog);
  const model = useStore((s) => s.model);
  // 从杆件打开 → 默认单元荷载; 从节点打开 → 默认节点集中力
  const [type, setType] = useState(
    dialog && dialog.element >= 0 ? "lateralUniformPressure" : "nodalForce"
  );
  const [direction, setDirection] = useState("y");
  const [value, setValue] = useState("-10000");
  const [element, setElement] = useState(dialog && dialog.element >= 0 ? String(dialog.element) : "0");
  const [position, setPosition] = useState("1");
  const [T0, setT0] = useState("0");
  const [T1, setT1] = useState("0");

  if (!dialog) return null;
  const st = useStore.getState();

  const isNodeLoad = type === "nodalForce" || type === "momentOnPoint" || type === "supportMove";
  const isElementLoad = !isNodeLoad; // 单元荷载(含温度)都需要选单元
  const isTemp = type === "temperature";

  function submit() {
    const v = Number(value);
    if (Number.isNaN(v)) return;
    // 节点荷载但当前没有有效节点(从杆件打开) → 提示切换
    if (isNodeLoad && dialog.node < 0) {
      st.setToast(t("toast.needNodeLoad"), true);
      return;
    }
    // 单元荷载但当前没有有效单元(从节点打开) → 提示切换
    if (isElementLoad && dialog.element < 0) {
      st.setToast(t("toast.needElemLoad"), true);
      return;
    }
    const params = { type, direction, value: v };
    if (isNodeLoad) {
      params.node = dialog.node;
    } else if (isElementLoad) {
      const eid = Number(element);
      if (Number.isNaN(eid) || !model.elements.some((e) => e.id === eid)) {
        st.setToast(t("toast.invalidElem"), true);
        return;
      }
      params.element = eid;
      const pos = Number(position);
      if (!Number.isNaN(pos)) params.position = pos;
      // 杆上集中力需要节点 = 单元起点 (后端按 iLoadedNode 定位)
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
        st.setToast(t("toast.tempNeedsElem"), true);
        return;
      }
      params.element = eid;
    }
    st.addLoad(params);
  }

  return (
    <div className="load-dialog">
      <div className="load-dialog-title">
        {t("loadDialog.title")} {isNodeLoad ? t("loadDialog.atNode", { n: dialog.node }) : ""}
      </div>
      <label className="field">
        <span>{t("loadDialog.type")}</span>
        <select value={type} onChange={(e) => setType(e.target.value)}>
          {Object.entries(LOAD_TYPES).map(([k]) => (
            <option key={k} value={k}>
              {t("loadType." + k)}
            </option>
          ))}
        </select>
      </label>
      {!isTemp && (
        <label className="field">
          <span>{t("loadDialog.direction")}</span>
          <select value={direction} onChange={(e) => setDirection(e.target.value)}>
            <option value="x">{t("loadDialog.dirX")}</option>
            <option value="y">{t("loadDialog.dirY")}</option>
            <option value="rz">{t("loadDialog.dirR")}</option>
          </select>
        </label>
      )}
      {isElementLoad && (
        <label className="field">
          <span>{t("loadDialog.element")}</span>
          <select value={element} onChange={(e) => setElement(e.target.value)}>
            {model.elements.map((e) => (
              <option key={e.id} value={e.id}>
                {t("loadDialog.elemOption", { id: e.id, a: e.nodeI, b: e.nodeJ })}
              </option>
            ))}
          </select>
        </label>
      )}
      {isTemp ? (
        <>
          <label className="field">
            <span>{t("loadDialog.t0")}</span>
            <input type="number" value={T0} onChange={(e) => setT0(e.target.value)} />
          </label>
          <label className="field">
            <span>{t("loadDialog.t1")}</span>
            <input type="number" value={T1} onChange={(e) => setT1(e.target.value)} />
          </label>
        </>
      ) : (
        <>
          <label className="field">
            <span>
              {type === "supportMove"
                ? t("loadDialog.valueDisp")
                : t("loadDialog.value") + " " + (isElementLoad ? t("loadDialog.valuePerM") : t("loadDialog.valueUnit"))}
            </span>
            <input type="number" value={value} onChange={(e) => setValue(e.target.value)} />
          </label>
          {isElementLoad && (
            <label className="field">
              <span>{t("loadDialog.position")}</span>
              <input type="number" value={position} onChange={(e) => setPosition(e.target.value)} />
            </label>
          )}
        </>
      )}
      <div className="load-dialog-actions">
        <button className="btn primary small" onClick={submit}>
          {t("loadDialog.confirm")}
        </button>
        <button className="btn small" onClick={() => st.setLoadDialog(null)}>
          {t("loadDialog.cancel")}
        </button>
      </div>
    </div>
  );
}
