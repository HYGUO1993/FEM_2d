#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LLM 案例接入脚本: 校验 + 求解 + 摘要。
用法:
  python scripts/llm_run.py model.json [result.json]
若 result.json 省略,输出到 model_result.json。
"""
import sys, os, json, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEMCLI = os.path.join(ROOT, "build", "bin", "femcli.exe")


def main():
    if len(sys.argv) < 2:
        print("用法: python scripts/llm_run.py <model.json> [result.json]")
        sys.exit(1)
    model_path = sys.argv[1]
    result_path = sys.argv[2] if len(sys.argv) > 2 else \
        os.path.splitext(model_path)[0] + "_result.json"

    # 1. 基本 JSON 校验
    try:
        with open(model_path, encoding="utf-8") as f:
            model = json.load(f)
    except Exception as e:
        print(f"[错误] 无法读取模型: {e}")
        sys.exit(1)

    # 2. femcli validate
    r = subprocess.run([FEMCLI, "validate", model_path],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[错误] 模型校验失败:\n{r.stdout}{r.stderr}")
        sys.exit(1)
    print(f"[OK] 校验: {r.stdout.strip()}")

    # 3. femcli solve
    r = subprocess.run([FEMCLI, "solve", model_path, "-o", result_path],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[错误] 求解失败:\n{r.stdout}{r.stderr}")
        sys.exit(1)

    with open(result_path, encoding="utf-8") as f:
        result = json.load(f)

    if result.get("status") != "ok":
        print(f"[错误] 求解器报错: {result.get('message')}")
        sys.exit(1)

    # 4. 摘要
    stats = result.get("stats", {})
    print(f"[OK] 求解完成: {stats.get('nodeCount')} 节点, "
          f"{stats.get('elementCount')} 单元, {stats.get('freeDOF')}/{stats.get('totalDOF')} DOF")
    print(f"[OK] 结果写入: {result_path}")
    print("\n位移摘要(非零):")
    for d in result.get("displacements", []):
        vals = {k: v for k, v in d.items() if k != "node" and abs(v) > 1e-12}
        if vals:
            print(f"  节点 {d['node']}: {vals}")
    print("\n反力:")
    for r_ in result.get("reactions", []):
        vals = {k: v for k, v in r_.items() if k != "node" and abs(v) > 1e-12}
        if vals:
            print(f"  节点 {r_['node']}: {vals}")


if __name__ == "__main__":
    main()
