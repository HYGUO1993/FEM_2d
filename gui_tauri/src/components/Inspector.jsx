import { useStore } from "../store.js";

function fmt(v) {
  if (typeof v !== "number" || Number.isNaN(v)) return String(v);
  const abs = Math.abs(v);
  if (abs >= 1e6 || (abs > 0 && abs < 1e-4)) return v.toExponential(3);
  return Number(v.toFixed(6)).toString();
}

export default function Inspector() {
  const selection = useStore((s) => s.selection);
  if (!selection) return <ModelPanel />;

  if (selection.type === "node") return <NodePanel id={selection.id} />;
  if (selection.type === "element") return <ElementPanel id={selection.id} />;
  if (selection.type === "load") return <LoadPanel id={selection.id} />;
  if (selection.type === "constraint") return <ConstraintPanel node={selection.id} />;
  return null;
}

// —— 节点属性 ——
function NodePanel({ id }) {
  const node = useStore((s) => s.model.nodes.find((n) => n.id === id));
  const updateNode = useStore((s) => s.updateNode);
  if (!node) return null;

  return (
    <div className="pane inspector">
      <h3>节点 #{id}</h3>
      <label className="field">
        <span>x (m)</span>
        <input
          type="number"
          step="0.1"
          value={node.x}
          onChange={(e) => {
            const v = parseFloat(e.target.value);
            if (!Number.isNaN(v)) updateNode(id, { x: v });
          }}
        />
      </label>
      <label className="field">
        <span>y (m)</span>
        <input
          type="number"
          step="0.1"
          value={node.y}
          onChange={(e) => {
            const v = parseFloat(e.target.value);
            if (!Number.isNaN(v)) updateNode(id, { y: v });
          }}
        />
      </label>
      <label className="field">
        <span>类型</span>
        <select value={node.type} onChange={(e) => updateNode(id, { type: e.target.value })}>
          <option value="frame">frame (刚架)</option>
          <option value="truss">truss (桁架)</option>
        </select>
      </label>
    </div>
  );
}

// —— 单元属性 ——
function ElementPanel({ id }) {
  const el = useStore((s) => s.model.elements.find((e) => e.id === id));
  const model = useStore((s) => s.model);
  const updateElement = useStore((s) => s.updateElement);
  if (!el) return null;

  return (
    <div className="pane inspector">
      <h3>杆件 #{id}</h3>
      <div className="readonly-row">
        <span>nodeI</span>
        <code>{el.nodeI}</code>
      </div>
      <div className="readonly-row">
        <span>nodeJ</span>
        <code>{el.nodeJ}</code>
      </div>
      <label className="field">
        <span>类型</span>
        <select value={el.type} onChange={(e) => updateElement(id, { type: e.target.value })}>
          <option value="frame">frame (刚架)</option>
          <option value="truss">truss (桁架)</option>
        </select>
      </label>
      <label className="field">
        <span>截面</span>
        <select
          value={el.section}
          onChange={(e) => updateElement(id, { section: Number(e.target.value) })}
        >
          {model.sections.map((sc) => (
            <option key={sc.id} value={sc.id}>
              #{sc.id} (A={sc.A}, Iz={sc.Iz})
            </option>
          ))}
        </select>
      </label>
      <label className="field">
        <span>材料</span>
        <select
          value={el.material}
          onChange={(e) => updateElement(id, { material: Number(e.target.value) })}
        >
          {model.materials.map((mt) => (
            <option key={mt.id} value={mt.id}>
              #{mt.id} (E={mt.E})
            </option>
          ))}
        </select>
      </label>
    </div>
  );
}

// —— 荷载属性 ——
function LoadPanel({ id }) {
  const ld = useStore((s) => s.model.loads.find((l) => l.id === id));
  const model = useStore((s) => s.model);
  const deleteLoad = useStore((s) => s.deleteLoad);
  if (!ld) return null;

  const byId = new Map(model.nodes.map((n) => [n.id, n]));
  const TYPE_NAMES = {
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

  return (
    <div className="pane inspector">
      <h3>荷载 #{id}</h3>
      <div className="readonly-row">
        <span>类型</span>
        <code>{TYPE_NAMES[ld.type] || ld.type}</code>
      </div>
      {ld.node >= 0 && (
        <div className="readonly-row">
          <span>节点</span>
          <code>
            {ld.node} ({byId.get(ld.node)?.x}, {byId.get(ld.node)?.y})
          </code>
        </div>
      )}
      {ld.element >= 0 && (
        <div className="readonly-row">
          <span>单元</span>
          <code>#{ld.element}</code>
        </div>
      )}
      {ld.type === "temperature" ? (
        <>
          <div className="readonly-row">
            <span>下表面 ΔT</span>
            <code>{ld.T0 ?? 0} ℃</code>
          </div>
          <div className="readonly-row">
            <span>上表面 ΔT</span>
            <code>{ld.T1 ?? 0} ℃</code>
          </div>
        </>
      ) : (
        <>
          <div className="readonly-row">
            <span>方向</span>
            <code>{ld.direction}</code>
          </div>
          <div className="readonly-row">
            <span>{ld.type === "supportMove" ? "位移值" : "数值"}</span>
            <code>{fmt(ld.value)} {ld.type === "supportMove" ? "m" : ""}</code>
          </div>
          {ld.position > 0 && (
            <div className="readonly-row">
              <span>位置/长度</span>
              <code>{ld.position} m</code>
            </div>
          )}
        </>
      )}
      <button className="btn danger small" onClick={() => deleteLoad(id)}>
        删除荷载
      </button>
    </div>
  );
}

// —— 约束属性 ——
function ConstraintPanel({ node }) {
  const model = useStore((s) => s.model);
  const toggleConstraint = useStore((s) => s.toggleConstraint);
  const constraint = model.constraints.find((c) => c.node === node);
  const dofs = constraint ? constraint.dofs : [];

  const opts = [
    { id: "ux", label: "ux (水平)" },
    { id: "uy", label: "uy (竖向)" },
    { id: "rz", label: "rz (转动)" },
  ];

  return (
    <div className="pane inspector">
      <h3>约束 @节点 {node}</h3>
      {opts.map((o) => (
        <label key={o.id} className="check-row">
          <input
            type="checkbox"
            checked={dofs.includes(o.id)}
            onChange={() => toggleConstraint(node, o.id)}
          />
          <span>{o.label}</span>
        </label>
      ))}
      <p className="hint">勾选 = 该自由度被约束 (位移 0)</p>
    </div>
  );
}

// —— 无选中: 模型全局设置 ——
function ModelPanel() {
  const model = useStore((s) => s.model);
  const updateMaterial = useStore((s) => s.updateMaterial);
  const updateSection = useStore((s) => s.updateSection);

  return (
    <div className="pane inspector">
      <h3>模型</h3>
      <div className="readonly-row">
        <span>标题</span>
        <code>{model.title}</code>
      </div>
      <div className="readonly-row">
        <span>节点 / 单元</span>
        <code>
          {model.nodes.length} / {model.elements.length}
        </code>
      </div>

      <h4 className="sub-heading">材料</h4>
      {model.materials.map((mt) => (
        <div key={mt.id} className="mat-row">
          <span className="mat-id">#{mt.id}</span>
          <label className="field">
            <span>E (Pa)</span>
            <input
              type="number"
              value={mt.E}
              onChange={(e) => {
                const v = parseFloat(e.target.value);
                if (!Number.isNaN(v)) updateMaterial(mt.id, { E: v });
              }}
            />
          </label>
          <div className="mat-inline">
            <label className="field">
              <span>μ</span>
              <input
                type="number"
                value={mt.mu}
                onChange={(e) => {
                  const v = parseFloat(e.target.value);
                  if (!Number.isNaN(v)) updateMaterial(mt.id, { mu: v });
                }}
              />
            </label>
            <label className="field">
              <span>α</span>
              <input
                type="number"
                value={mt.alpha}
                onChange={(e) => {
                  const v = parseFloat(e.target.value);
                  if (!Number.isNaN(v)) updateMaterial(mt.id, { alpha: v });
                }}
              />
            </label>
          </div>
        </div>
      ))}

      <h4 className="sub-heading">截面</h4>
      {model.sections.map((sc) => (
        <div key={sc.id} className="mat-row">
          <span className="mat-id">#{sc.id}</span>
          <label className="field">
            <span>A (m²)</span>
            <input
              type="number"
              value={sc.A}
              onChange={(e) => {
                const v = parseFloat(e.target.value);
                if (!Number.isNaN(v)) updateSection(sc.id, { A: v });
              }}
            />
          </label>
          <div className="mat-inline">
            <label className="field">
              <span>Iz (m⁴)</span>
              <input
                type="number"
                value={sc.Iz}
                onChange={(e) => {
                  const v = parseFloat(e.target.value);
                  if (!Number.isNaN(v)) updateSection(sc.id, { Iz: v });
                }}
              />
            </label>
            <label className="field">
              <span>h (m)</span>
              <input
                type="number"
                value={sc.height}
                onChange={(e) => {
                  const v = parseFloat(e.target.value);
                  if (!Number.isNaN(v)) updateSection(sc.id, { height: v });
                }}
              />
            </label>
          </div>
        </div>
      ))}
    </div>
  );
}
