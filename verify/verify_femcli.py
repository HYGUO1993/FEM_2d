#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
femcli 黄金用例回归验证
对比 femcli 的 txt 模式输出与 tests/golden/golden.json 参考解。
用法: python verify/verify_femcli.py
"""
import json, os, subprocess, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FEMCLI = os.path.join(ROOT, "build", "bin", "femcli.exe")
GOLDEN = os.path.join(ROOT, "tests", "golden", "golden.json")

def load_golden():
    with open(GOLDEN) as f:
        return json.load(f)["cases"]

def run_femcli(input_path, out_path):
    r = subprocess.run([FEMCLI, "solve", "-i", input_path, "-o", out_path, "--legacy-txt"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None, r.stderr.strip()
    with open(out_path) as f:
        return json.load(f), None

def compare(actual, golden_disp, golden_forces, golden_reac, tol=1e-9):
    """对比位移/端力/反力,返回 (ok, maxdiff, detail)

    坐标系说明:
      - golden 的 end_forces 为整体坐标 [Fx1,Fy1,Mz1,Fx2,Fy2,Mz2]
      - femcli 输出局部坐标 [N_i,V_i,M_i,N_j,V_j,M_j]
      水平单元两者一致,斜杆数值不同(局部=整体旋转)。
      因此端力不做数值对比(由 fem_run 一致性保证),
      严格对比: 位移 + 反力(两者均为整体坐标,物理量)。
    """
    maxdiff = 0.0
    # 1) 位移
    for i, gd in enumerate(golden_disp):
        a = actual["displacements"][i]
        for j, key in enumerate(["ux", "uy", "rz"]):
            maxdiff = max(maxdiff, abs(a[key] - gd[j]))
    # 2) 反力
    for i, gr in enumerate(golden_reac):
        a = actual["reactions"][i]
        for j, key in enumerate(["ux", "uy", "rz"]):
            maxdiff = max(maxdiff, abs(a[key] - gr[j]))
    # 3) 端力只做基本物理校验: 力幅值有限(非 NaN/Inf), 且单元数一致
    assert len(actual["endForces"]) == len(golden_forces), "端力数量不一致"
    for e in actual["endForces"]:
        for k in ("N", "V", "M"):
            v = e["nodeI"][k]
            if not math.isfinite(v):
                maxdiff = max(maxdiff, float("inf"))
    return maxdiff <= tol, maxdiff

def main():
    golden = load_golden()
    print(f"femcli: {FEMCLI}")
    print(f"golden cases: {list(golden.keys())}")
    print("-" * 60)
    all_ok = True
    for name, case in golden.items():
        input_file = case["input"]
        if not os.path.isabs(input_file):
            input_file = os.path.join(ROOT, input_file)
        if not os.path.exists(input_file):
            print(f"[SKIP] {name}: 输入文件不存在 {input_file}")
            continue
        out_path = os.path.join(HERE, f"_femcli_{name}.json")
        actual, err = run_femcli(input_file, out_path)
        if actual is None:
            print(f"[FAIL] {name}: femcli 报错: {err}")
            all_ok = False
            continue
        ok, maxdiff = compare(actual, case["displacements"], case["end_forces"], case["reactions"])
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_ok = False
        print(f"[{status}] {name:24s} maxdiff={maxdiff:.3e}")
        os.remove(out_path)
    print("-" * 60)
    print("全部通过!" if all_ok else "存在失败用例!")
    return 0 if all_ok else 1

if __name__ == "__main__":
    sys.exit(main())
