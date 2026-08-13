import { useStore } from "../store.js";
import { t } from "../i18n.js";

function fmt(v) {
  if (typeof v !== "number" || Number.isNaN(v)) return String(v);
  if (v === 0) return "0";
  const abs = Math.abs(v);
  if (abs >= 1e6 || (abs > 0 && abs < 1e-4)) return v.toExponential(3);
  return Number(v.toFixed(6)).toString();
}

export default function ResultsPanel() {
  const results = useStore((s) => s.results);
  const solved = useStore((s) => s.solved);
  const solveTime = useStore((s) => s.solveTime);
  const solving = useStore((s) => s.solving);
  const deformScale = useStore((s) => s.deformScale);
  const setDeformScale = useStore((s) => s.setDeformScale);
  const solve = useStore((s) => s.solve);

  return (
    <div className="pane results-panel">
      <div className="pane-head">
        <h3>{t("solveResults")}</h3>
        {solveTime && <span className="muted">{solveTime}s</span>}
      </div>
      <button className="btn primary block" onClick={() => solve()} disabled={solving}>
        {solving ? t("results.solving") : t("results.solveBtn")}
      </button>

      {!results && <p className="placeholder">{t("results.placeholder")}</p>}

      {results && results.status === "error" && (
        <div className="error-box">
          <h4>{t("results.failed")}</h4>
          <p>{results.message || t("results.unknownErr")}</p>
        </div>
      )}

      {results && results.status === "ok" && (
        <>
          {/* stats 卡片 */}
          <div className="stats-row">
            {[
              ["results.statNodes", results.stats?.nodeCount],
              ["results.statElems", results.stats?.elementCount],
              ["results.statDof", results.stats?.freeDOF != null ? `${results.stats.freeDOF}/${results.stats.totalDOF}` : "-"],
            ].map(([label, val]) => (
              <div className="stat-card" key={label}>
                <div className="label">{t(label)}</div>
                <div className="value">{val ?? "-"}</div>
              </div>
            ))}
          </div>

          {/* 内力极值汇总 */}
          {(results.endForces || []).length > 0 && (
            <ExtremeCard results={results} />
          )}

          {/* 变形图倍率 */}
          {solved && (
            <div className="deform-row">
              <span>{t("results.deformScale")}</span>
              <input
                type="range"
                min="0"
                max="500"
                step="1"
                value={deformScale}
                onChange={(e) => setDeformScale(Number(e.target.value))}
              />
              <code>{Math.round(deformScale)}x</code>
            </div>
          )}

          {/* 位移 */}
          <h4 className="sub-heading">{t("results.displacements")}</h4>
          <table className="res-table">
            <thead>
              <tr>
                <th>{t("results.thNode")}</th>
                <th>ux</th>
                <th>uy</th>
                <th>rz</th>
              </tr>
            </thead>
            <tbody>
              {(results.displacements || []).slice(0, 500).map((d) => (
                <tr key={d.node}>
                  <td>{d.node}</td>
                  <td className={`num ${d.ux === 0 ? "zero" : ""}`}>{fmt(d.ux)}</td>
                  <td className={`num ${d.uy === 0 ? "zero" : ""}`}>{fmt(d.uy)}</td>
                  <td className={`num ${d.rz === 0 ? "zero" : ""}`}>{fmt(d.rz)}</td>
                </tr>
              ))}
            </tbody>
          </table>

          {/* 端力 */}
          <h4 className="sub-heading">{t("results.endForces")}</h4>
          <p className="hint">{t("results.momentNote")}</p>
          <table className="res-table">
            <thead>
              <tr>
                <th>{t("results.thElement")}</th>
                <th>N_i</th>
                <th>V_i</th>
                <th>M_i</th>
                <th>N_j</th>
                <th>V_j</th>
                <th>M_j</th>
              </tr>
            </thead>
            <tbody>
              {(results.endForces || []).slice(0, 500).map((e) => (
                <tr key={e.element}>
                  <td>{e.element}</td>
                  <td className="num">{fmt(e.nodeI.N)}</td>
                  <td className="num">{fmt(e.nodeI.V)}</td>
                  <td className="num">{fmt(e.nodeI.M)}</td>
                  <td className="num">{fmt(e.nodeJ.N)}</td>
                  <td className="num">{fmt(e.nodeJ.V)}</td>
                  <td className="num">{fmt(e.nodeJ.M)}</td>
                </tr>
              ))}
            </tbody>
          </table>

          {/* 反力 */}
          <h4 className="sub-heading">{t("results.reactions")}</h4>
          <table className="res-table">
            <thead>
              <tr>
                <th>{t("results.thNode")}</th>
                <th>Rx</th>
                <th>Ry</th>
                <th>Rz</th>
              </tr>
            </thead>
            <tbody>
              {(results.reactions || []).slice(0, 500).map((r) => (
                <tr key={r.node}>
                  <td>{r.node}</td>
                  <td className={`num ${r.ux === 0 ? "zero" : ""}`}>{fmt(r.ux)}</td>
                  <td className={`num ${r.uy === 0 ? "zero" : ""}`}>{fmt(r.uy)}</td>
                  <td className={`num ${r.rz === 0 ? "zero" : ""}`}>{fmt(r.rz)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </div>
  );
}

/** 内力极值汇总: 扫描所有单元端力, 输出 N/V/M 的 max|值| 及所在单元 */
function ExtremeCard({ results }) {
  const ef = results.endForces || [];
  if (!ef.length) return null;
  const keys = ["N", "V", "M"];
  const stats = keys.map((k) => {
    let maxAbs = 0;
    let maxAbsVal = 0;
    let elem = "-";
    for (const e of ef) {
      for (const side of ["nodeI", "nodeJ"]) {
        const v = e[side]?.[k];
        if (typeof v === "number" && Math.abs(v) > maxAbs) {
          maxAbs = Math.abs(v);
          maxAbsVal = v;
          elem = e.element;
        }
      }
    }
    return { k, maxAbsVal, elem };
  });

  return (
    <div className="extreme-row">
      {stats.map((s) => (
        <div className="stat-card" key={s.k}>
          <div className="label">{t("results.maxLabel", { k: s.k, e: s.elem })}</div>
          <div className={`value ${s.maxAbsVal < 0 ? "neg" : ""}`}>{fmt(s.maxAbsVal)}</div>
        </div>
      ))}
    </div>
  );
}
