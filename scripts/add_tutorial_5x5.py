#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 examples/frame_5x5.json 转成 tutorialProjects.js 内联条目并插入"""
import json, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with open(os.path.join(ROOT, "examples", "frame_5x5.json"), encoding="utf-8") as f:
    model = json.load(f)
compact = json.dumps(model, ensure_ascii=False, separators=(",", ":"))
entry = '  { id: "frame_5x5", desc: "5x5 刚框架 · 水平+竖向 10kN (5跨x5层)", json: "%s" },\n' % compact.replace('"', '\\"')

path = os.path.join(ROOT, "gui_tauri", "src", "tutorialProjects.js")
src = open(path, encoding="utf-8").read()
marker = '  { id: "cantilever_point"'
assert marker in src, "marker not found"
src = src.replace(marker, entry + marker)
open(path, "w", encoding="utf-8").write(src)
print("Inserted tutorial entry, len:", len(entry))
