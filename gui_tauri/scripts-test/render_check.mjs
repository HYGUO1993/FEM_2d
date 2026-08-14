// render.js 单元测试: 内力分布 + Hermite 变形图几何验证 (Node 直接运行)
// 用法: node gui_tauri/scripts-test/render_check.mjs
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
const require = createRequire(import.meta.url);
import { computeForceDiagram, computeDeformed, loadContribAt } from "../src/canvas/render.js";

const ROOT = new URL("../..", import.meta.url).pathname.replace(/^\//, "").replace(/\//g, "\\");

function load(path) {
  return JSON.parse(readFileSync(ROOT + path, "utf-8"));
}

// 解析 path d 字符串 → 世界坐标数组 (view=null: sx=x, sy=-y)
function decode(d) {
  return d
    .replace(/^M /, "")
    .split(" L ")
    .map((p) => {
      const [sx, sy] = p.split(",").map(Number);
      return { x: sx, y: -sy };
    });
}

let fails = 0;
function check(name, cond, detail) {
  if (!cond) fails++;
  console.log(`[${cond ? "PASS" : "FAIL"}] ${name}${detail ? "  " + detail : ""}`);
}
function approx(a, b, tol = 1e-6) {
  return Math.abs(a - b) <= tol;
}

// ============ 1) loadContribAt: 荷载贡献公式 ============
{
  const udl = [{ type: "lateralUniformPressure", value: -20000, position: 6 }];
  const c0 = loadContribAt(udl, 6, 0);
  const c3 = loadContribAt(udl, 6, 3);
  const c6 = loadContribAt(udl, 6, 6);
  check("UDL V(0)=0", approx(c0.V, 0));
  check("UDL V(3)=-60000", approx(c3.V, -60000));
  check("UDL V(6)=-120000", approx(c6.V, -120000));
  check("UDL M(3)=-90000", approx(c3.M, -90000), `M=${c3.M}`);
  check("UDL M(6)=-360000", approx(c6.M, -360000));
  const pt = [{ type: "lateralForce", value: -3000, position: 2 }];
  const p0 = loadContribAt(pt, 6, 1);
  const p1 = loadContribAt(pt, 6, 2);
  const p2 = loadContribAt(pt, 6, 3);
  check("P V jump", approx(p0.V, 0) && approx(p1.V, -3000) && approx(p2.V, -3000));
  check("P M after", approx(p2.M, -3000), `M=${p2.M}`);
  const lin = [{ type: "lateralLinearlyPressure", value: -1000, position: 4 }];
  const l4 = loadContribAt(lin, 4, 4);
  check("linear V(4)=-2000", approx(l4.V, -2000), `V=${l4.V}`);
  check("linear M(4)=-2666.67", approx(l4.M, -2666.6667, 1e-3), `M=${l4.M}`);
  const ax = [{ type: "axialPressure", value: 20000, position: 4 }];
  check("axial N(4)=80000", approx(loadContribAt(ax, 4, 4).N, 80000));
  const axf = [{ type: "axialForce", value: 100000, position: 4 }];
  check("axialF N(4)=100000", approx(loadContribAt(axf, 4, 4).N, 100000));
}

// ============ 2) 门式刚架 UDL: M 图 (反弯点!) / V 图 (线性) ============
{
  const model = load("verify/portal_udl_test.json");
  const { execSync } = require("node:child_process");
  execSync(
    `"${ROOT}build\\bin\\Release\\femcli.exe" solve "${ROOT}verify\\portal_udl_test.json" -o "${ROOT}verify\\_render_tmp.json"`,
    { encoding: "utf-8" }
  );
  const results = JSON.parse(readFileSync(ROOT + "verify\\_render_tmp.json", "utf-8"));
  const M_i = 41466.69226339029;
  const V_i = 60000;
  const q = -20000;
  const L = 6;
  const Mmid = M_i - V_i * 3 - (q * 9) / 2; // 理论跨中 M
  const disc = Math.sqrt(V_i * V_i + 2 * q * M_i);
  const inf0 = (V_i - disc) / -q; // M(x)=0 两根
  const inf1 = (V_i + disc) / -q;
  const at = (samples, x) => samples.find((p) => approx(p.x, x, 1e-9))?.v;
  const near = (samples, x, tol) => {
    const hit = samples.filter((p) => Math.abs(p.x - x) < tol);
    return hit.length ? hit[0].v : NaN;
  };

  const mPath = computeForceDiagram(model, results, "M", null).find((p) => p.key === "fg-M-2");
  const vPath = computeForceDiagram(model, results, "V", null).find((p) => p.key === "fg-V-2");
  check("M/V 图存在", !!mPath && !!vPath);
  if (mPath && vPath) {
    const sM = mPath.samples;
    const sV = vPath.samples;
    // 端值
    check("M 两端同号且幅值相等", approx(sM[0].v, M_i, 1e-6) && approx(sM[sM.length - 1].v, M_i, 1e-6),
      `M_end=${sM[0].v.toFixed(1)}`);
    // 反弯点: M=0 的两根, 用线性插值找采样间的零点
    const crossings = [];
    for (let i = 1; i < sM.length; i++) {
      const v0 = sM[i - 1].v;
      const v1 = sM[i].v;
      if (v0 * v1 < 0) crossings.push(sM[i - 1].x + ((sM[i].x - sM[i - 1].x) * v0) / (v0 - v1));
    }
    check("M 图含 2 个反弯点", crossings.length === 2, `crossings=${crossings.map((x) => x.toFixed(3)).join(",")} (理论 ${inf0.toFixed(3)},${inf1.toFixed(3)})`);
    if (crossings.length === 2) {
      check("反弯点 x 位置准确", approx(crossings[0], inf0, 0.02) && approx(crossings[1], inf1, 0.02),
        `x=${crossings.map((x) => x.toFixed(3)).join(",")}`);
    }
    // 跨中弯矩
    const vMid = near(sM, 3, 0.2);
    check("M 跨中值匹配理论", approx(vMid, Mmid, Math.abs(Mmid) * 1e-6), `Mmid=${vMid.toFixed(1)} (理论 ${Mmid.toFixed(1)})`);
    // V: 线性 +60000 → -60000, 跨中 0
    check("V 端值 +60000/-60000", approx(sV[0].v, V_i, 1e-6) && approx(sV[sV.length - 1].v, -V_i, 1e-6),
      `V=${sV[0].v.toFixed(1)}→${sV[sV.length - 1].v.toFixed(1)}`);
    const vMidV = near(sV, 3, 0.2);
    check("V 跨中为 0", approx(vMidV, 0, 1e-6), `v_mid=${vMidV.toFixed(2)}`);
    // 几何: path 必须沿杆轴展开 (x 覆盖单元起点→终点), 防止缩成起点处竖线
    const beamPts = decode(mPath.d);
    const xs = beamPts.map((p) => p.x);
    const xMin = Math.min(...xs);
    const xMax = Math.max(...xs);
    check("M 图沿杆轴展开 (x∈[0,6])", approx(xMin, 0, 1e-6) && approx(xMax, 6, 1e-6),
      `x=[${xMin.toFixed(2)},${xMax.toFixed(2)}] (期望 [0,6])`);
    // 柱: 无荷载, M 线性 0 → 顶端 41466.69
    const colM = computeForceDiagram(model, results, "M", null).find((p) => p.key === "fg-M-0");
    const sC = colM.samples;
    check("柱 M 底部≈0 顶端=41466.69",
      approx(sC[0].v, 0, 1e-6) && approx(sC[sC.length - 1].v, 41466.69, 1e-2),
      `v=${sC[sC.length - 1].v.toFixed(1)}`);
  }
}

// ============ 3) 简支梁 UDL: M 抛物线跨中 -5000, V 三角 +10000→-10000 ============
{
  const model = load("examples/simply_supported_udl.json");
  const results = JSON.parse(readFileSync(ROOT + "gui_tauri\\scripts-test\\_render_ssb.json", "utf-8"));
  const at = (samples, x) => samples.find((p) => approx(p.x, x, 1e-9))?.v;
  const m0 = computeForceDiagram(model, results, "M", null).find((p) => p.key === "fg-M-0").samples;
  const m1 = computeForceDiagram(model, results, "M", null).find((p) => p.key === "fg-M-1").samples;
  const v0 = computeForceDiagram(model, results, "V", null).find((p) => p.key === "fg-V-0").samples;
  const v1 = computeForceDiagram(model, results, "V", null).find((p) => p.key === "fg-V-1").samples;
  // elem0: 铰支端 M=0 → j 端 -3750 (跨节点连续值)
  check("SSB M 铰支端=0", approx(at(m0, 0), 0, 1e-6), `M=${at(m0, 0)}`);
  check("SSB M 节点连续 (elem0 j = elem1 i = -3750)",
    approx(at(m0, 0.5), -3750, 1e-6) && approx(at(m1, 0), -3750, 1e-6),
    `M=${at(m0, 0.5).toFixed(1)}/${at(m1, 0).toFixed(1)}`);
  // elem1: 跨中 (j 端, 局部 x=0.5) M = -5000 (理论峰值)
  check("SSB M 跨中=-5000", approx(at(m1, 0.5), -5000, 1e-6), `Mmid=${at(m1, 0.5).toFixed(1)}`);
  // V: elem0 +10000 → +5000 (跨节点连续), 跨中 0
  check("SSB V 端值 +10000", approx(at(v0, 0), 10000, 1e-3), `V=${at(v0, 0).toFixed(0)}`);
  check("SSB V 四分之一跨 =+5000", approx(at(v0, 0.5), 5000, 1e-3), `V=${at(v0, 0.5).toFixed(0)}`);
  check("SSB V 跨中=0", approx(at(v1, 0.5), 0, 1e-6), `V=${at(v1, 0.5).toFixed(2)}`);
}

// ============ 4) Hermite 变形图: 悬臂梁端弯矩 (抛物线形) ============
{
  const model = load("verify/moment_test.json");
  const results = JSON.parse(readFileSync(ROOT + "gui_tauri\\scripts-test\\_render_moment.json", "utf-8"));
  const defs = computeDeformed(model, results, 1, null);
  const d = defs.find((x) => x.key === "def-0");
  const pts = decode(d.d);
  const vEnd = pts[pts.length - 1].y - 0;
  const vMid = pts[Math.round(pts.length / 2) - 1].y - 0;
  check("变形图抛物线 (v_mid/v_end=0.25)", approx(vMid / vEnd, 0.25, 1e-9), `ratio=${(vMid / vEnd).toFixed(6)}`);
}

console.log(fails === 0 ? "\n全部通过!" : `\n${fails} 项失败!`);
process.exit(fails === 0 ? 0 : 1);
