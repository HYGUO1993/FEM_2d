#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 5x5 刚框架模型: 5跨 x 5层, 每跨/层 1m, 底部铰支, 顶层水平+竖向 10kN"""
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NX, NY = 6, 6  # 节点列数/层数 (5跨 x 5层)
SPAN = 1.0

nodes = []
for j in range(NY):
    for i in range(NX):
        nodes.append({"id": j * NX + i, "type": "frame", "x": i * SPAN, "y": j * SPAN})

constraints = [{"node": i, "dofs": ["ux", "uy"]} for i in range(NX)]  # 底层铰支

elements = []
eid = 0
# 水平杆: 每层每跨
for j in range(NY):
    for i in range(NX - 1):
        elements.append({"id": eid, "type": "frame", "nodeI": j * NX + i, "nodeJ": j * NX + i + 1, "section": 0, "material": 0})
        eid += 1
# 竖向杆: 每列每层
for i in range(NX):
    for j in range(NY - 1):
        elements.append({"id": eid, "type": "frame", "nodeI": j * NX + i, "nodeJ": (j + 1) * NX + i, "section": 0, "material": 0})
        eid += 1

loads = []
# 顶层水平荷载 10kN (+x) — 顶层 6 节点
top = [5 * NX + i for i in range(NX)]
for n in top:
    loads.append({"type": "nodalForce", "direction": "x", "value": 10000, "node": n})
# 顶层竖向荷载 10kN (向下) — 顶层中间 4 节点
for i in range(1, NX - 1):
    loads.append({"type": "nodalForce", "direction": "y", "value": -10000, "node": 5 * NX + i})

model = {
    "schemaVersion": "1.0",
    "title": "5x5 刚框架 - 水平+竖向 10kN (5跨x5层)",
    "solver": "builtin",
    "nodes": nodes,
    "constraints": constraints,
    "elements": elements,
    "materials": [{"id": 0, "E": 210000000000, "mu": 0.3, "alpha": 0}],
    "sections": [{"id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.2}],
    "loads": loads,
}

out = os.path.join(ROOT, "examples", "frame_5x5.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(model, f, ensure_ascii=False, indent=2)
print(f"Wrote {out}: {len(nodes)} 节点, {len(elements)} 单元, {len(loads)} 荷载")
