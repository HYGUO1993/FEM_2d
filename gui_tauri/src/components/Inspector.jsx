import { useStore } from "../store.js";
import { t } from "../i18n.js";

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
      <h3>{t("inspector.nodeTitle", { id })}</h3>
      <label className="field">
        <span>{t("inspector.x")}</span>
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
        <span>{t("inspector.y")}</span>
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
        <span>{t("inspector.type")}</span>
        <select value={node.type} onChange={(e) => updateNode(id, { type: e.target.value })}>
          <option value="frame">{t("inspector.typeFrame")}</option>
          <option value="truss">{t("inspector.typeTruss")}</option>
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
      <h3>{t("inspector.elemTitle", { id })}</h3>
      <div className="readonly-row">
        <span>nodeI</span>
        <code>{el.nodeI}</code>
      </div>
      <div className="readonly-row">
        <span>nodeJ</span>
        <code>{el.nodeJ}</code>
      </div>
      <label className="field">
        <span>{t("inspector.type")}</span>
        <select value={el.type} onChange={(e) => updateElement(id, { type: e.target.value })}>
          <option value="frame">{t("inspector.typeFrame")}</option>
          <option value="truss">{t("inspector.typeTruss")}</option>
        </select>
      </label>
      <label className="field">
        <span>{t("inspector.section")}</span>
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
        <span>{t("inspector.material")}</span>
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

  return (
    <div className="pane inspector">
      <h3>{t("inspector.loadTitle", { id })}</h3>
      <div className="readonly-row">
        <span>{t("inspector.type")}</span>
        <code>{t("loadType." + ld.type) !== "loadType." + ld.type ? t("loadType." + ld.type) : ld.type}</code>
      </div>
      {ld.node >= 0 && (
        <div className="readonly-row">
          <span>{t("inspector.node")}</span>
          <code>
            {ld.node} ({byId.get(ld.node)?.x}, {byId.get(ld.node)?.y})
          </code>
        </div>
      )}
      {ld.element >= 0 && (
        <div className="readonly-row">
          <span>{t("inspector.element")}</span>
          <code>#{ld.element}</code>
        </div>
      )}
      {ld.type === "temperature" ? (
        <>
          <div className="readonly-row">
            <span>{t("inspector.t0")}</span>
            <code>{ld.T0 ?? 0} ℃</code>
          </div>
          <div className="readonly-row">
            <span>{t("inspector.t1")}</span>
            <code>{ld.T1 ?? 0} ℃</code>
          </div>
        </>
      ) : (
        <>
          <div className="readonly-row">
            <span>{t("inspector.direction")}</span>
            <code>{ld.direction}</code>
          </div>
          <div className="readonly-row">
            <span>{ld.type === "supportMove" ? t("inspector.dispValue") : t("inspector.value")}</span>
            <code>{fmt(ld.value)} {ld.type === "supportMove" ? "m" : ""}</code>
          </div>
          {ld.position > 0 && (
            <div className="readonly-row">
              <span>{t("inspector.position")}</span>
              <code>{ld.position} m</code>
            </div>
          )}
        </>
      )}
      <button className="btn danger small" onClick={() => deleteLoad(id)}>
        {t("inspector.deleteLoad")}
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
    { id: "ux", label: t("inspector.ux") },
    { id: "uy", label: t("inspector.uy") },
    { id: "rz", label: t("inspector.rz") },
  ];

  return (
    <div className="pane inspector">
      <h3>{t("inspector.constraintTitle", { n: node })}</h3>
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
      <p className="hint">{t("inspector.constraintHint")}</p>
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
      <h3>{t("inspector.modelTitle")}</h3>
      <div className="readonly-row">
        <span>{t("inspector.title")}</span>
        <code>{model.title}</code>
      </div>
      <div className="readonly-row">
        <span>{t("inspector.nodesElems")}</span>
        <code>
          {model.nodes.length} / {model.elements.length}
        </code>
      </div>

      <h4 className="sub-heading">{t("inspector.material")}</h4>
      {model.materials.map((mt) => (
        <div key={mt.id} className="mat-row">
          <span className="mat-id">#{mt.id}</span>
          <label className="field">
            <span>{t("inspector.E")}</span>
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
              <span>{t("inspector.mu")}</span>
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
              <span>{t("inspector.alpha")}</span>
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

      <h4 className="sub-heading">{t("inspector.section")}</h4>
      {model.sections.map((sc) => (
        <div key={sc.id} className="mat-row">
          <span className="mat-id">#{sc.id}</span>
          <label className="field">
            <span>{t("inspector.A")}</span>
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
              <span>{t("inspector.Iz")}</span>
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
              <span>{t("inspector.h")}</span>
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
