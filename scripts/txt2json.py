#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
txt → JSON 转换工具
把 FEM_2d 传统 .txt 输入文件转成统一 JSON 模型(femcli 格式)。
用法: python scripts/txt2json.py test_beam.txt -o beam.json

注意: 按行解析(与 C++ 读取逻辑一致)——单元行材料索引可选,
      仅当该行有第 5 个 token 时才解析。
"""
import sys, os, json

# txt 数值载荷类型 → JSON 字符串
LOAD_TYPES = {
    1: "nodalForce",
    2: "lateralForce",
    3: "lateralUniformPressure",
    4: "momentOnPoint",
    5: "lateralLinearlyPressure",
    6: "axialPressure",
    7: "axialForce",
    8: "momentOnBeam",
    9: "temperature",
    10: "supportMove",
}
DIRS = {0: "x", 1: "y", 2: "rz"}
NODE_TYPES = {1: "truss", 2: "frame"}
ELEM_TYPES = {1: "truss", 2: "frame"}


def txt_to_json(txt_path):
    with open(txt_path, "r", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if ln.strip()]

    def toks_of(ln):
        return ln.split()

    # 总控行
    hdr = toks_of(lines[0])
    nNode, nCons, nElem, nMat, nSec, nLoad = map(int, hdr[0:6])
    ln = 1

    nodes = []
    for i in range(nNode):
        t, x, y = toks_of(lines[ln])
        nodes.append({"id": i, "type": NODE_TYPES.get(int(t), "frame"),
                      "x": float(x), "y": float(y)})
        ln += 1

    cons = []
    for i in range(nCons):
        n, dx, dy, dr = toks_of(lines[ln])
        dofs = []
        if int(dx) < 0: dofs.append("ux")
        if int(dy) < 0: dofs.append("uy")
        if int(dr) < 0: dofs.append("rz")
        cons.append({"node": int(n), "dofs": dofs})
        ln += 1

    elems = []
    for i in range(nElem):
        t = toks_of(lines[ln])
        mat = 0
        if len(t) >= 5:  # 材料索引可选
            mat = int(t[4])
        elems.append({"id": i, "type": ELEM_TYPES.get(int(t[0]), "frame"),
                      "nodeI": int(t[1]), "nodeJ": int(t[2]),
                      "section": int(t[3]), "material": mat})
        ln += 1

    mats = []
    for i in range(nMat):
        e, mu, alpha = toks_of(lines[ln])
        mats.append({"id": i, "E": float(e), "mu": float(mu), "alpha": float(alpha)})
        ln += 1

    sects = []
    for i in range(nSec):
        a, iz, h = toks_of(lines[ln])
        sects.append({"id": i, "A": float(a), "Iz": float(iz), "height": float(h)})
        ln += 1

    loads = []
    for i in range(nLoad):
        v = toks_of(lines[ln])
        lt = int(v[0])
        loads.append({
            "type": LOAD_TYPES.get(lt, str(lt)),
            "direction": DIRS.get(int(v[1]), "y"),
            "value": float(v[2]),
            "element": int(v[3]),
            "node": int(v[4]),
            "position": float(v[5]),
            "T0": float(v[6]) if len(v) > 6 else 0.0,
            "T1": float(v[7]) if len(v) > 7 else 0.0,
        })
        ln += 1

    return {
        "schemaVersion": "1.0",
        "title": os.path.splitext(os.path.basename(txt_path))[0],
        "solver": "builtin",
        "nodes": nodes,
        "constraints": cons,
        "elements": elems,
        "materials": mats,
        "sections": sects,
        "loads": loads,
    }


def main():
    if len(sys.argv) < 2:
        print("用法: python scripts/txt2json.py <input.txt> [-o output.json]")
        sys.exit(1)
    src = sys.argv[1]
    out = "model.json"
    if "-o" in sys.argv:
        out = sys.argv[sys.argv.index("-o") + 1]
    model = txt_to_json(src)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(model, f, ensure_ascii=False, indent=2)
    print(f"转换完成: {src} -> {out} ({len(model['nodes'])} 节点, "
          f"{len(model['elements'])} 单元, {len(model['loads'])} 荷载)")


if __name__ == "__main__":
    main()
