#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
femcli vs 独立参考求解器 (portal_udl_ref.solve) 对比验证。
用法: python verify/compare_femcli_ref.py <model.json> [<model2.json> ...]
输出: 位移/反力/杆端力逐项 maxdiff, 全 PASS 时退出码 0。
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FEMCLI = os.path.join(ROOT, "build", "bin", "Release", "femcli.exe")
sys.path.insert(0, HERE)
import portal_udl_ref as ref

def run_femcli(model_path):
    out = os.path.join(HERE, "_cmp_out.json")
    r = subprocess.run([FEMCLI, "solve", model_path, "-o", out],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None, r.stderr.strip()
    with open(out, encoding="utf-8") as f:
        return json.load(f), None
    os.remove(out)

def maxdiff(a, b):
    if a is None or b is None:
        return float("inf")
    d = 0.0
    for key in ("ux", "uy", "rz"):
        d = max(d, abs((a.get(key) or 0.0) - (b.get(key) or 0.0)))
    return d

def main():
    cases = sys.argv[1:] or [os.path.join(ROOT, "examples", "portal_frame.json")]
    all_ok = True
    for path in cases:
        with open(path, encoding="utf-8") as f:
            model = json.load(f)
        disp, endf, reac = ref.solve(model)
        actual, err = run_femcli(path)
        if actual is None:
            print("[FAIL] %s: femcli 报错: %s" % (os.path.basename(path), err))
            all_ok = False
            continue
        d = 0.0
        for r in disp:
            for a in actual["displacements"]:
                if a["node"] == r["node"]:
                    d = max(d, maxdiff(a, r)); break
        for r in reac:
            for a in actual["reactions"]:
                if a["node"] == r["node"]:
                    d = max(d, maxdiff(a, r)); break
        for e in endf:
            for a in actual["endForces"]:
                if a["element"] == e["element"]:
                    for side in ("nodeI", "nodeJ"):
                        for k in ("N", "V", "M"):
                            d = max(d, abs((a[side].get(k) or 0.0) - (e[side].get(k) or 0.0)))
                    break
        ok = d < 1e-6
        all_ok = all_ok and ok
        print("[%s] %-28s maxdiff=%.3e" % ("PASS" if ok else "FAIL", os.path.basename(path), d))
    print("全部通过!" if all_ok else "存在失败!")
    return 0 if all_ok else 1

if __name__ == "__main__":
    sys.exit(main())
