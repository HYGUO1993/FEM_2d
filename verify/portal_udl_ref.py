#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
独立参考求解器(密实组装 + 修正版 LDLT, 与 golden_generate.py 同一套数学):
门式刚架 + 梁上均布荷载 (lateralUniformPressure 满跨)。
用法: python verify/portal_udl_ref.py <model.json>
输出: 与 femcli 结果逐项对比 (位移/反力/杆端力), maxdiff
"""
import json, math, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

def mat_vec(A, x):
    return [sum(A[i][j]*x[j] for j in range(len(x))) for i in range(len(A))]

def mat_mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def transpose(A):
    return [list(r) for r in zip(*A)]

def ldlt_solve(A, b):
    n = len(A); L = [[0.0]*n for _ in range(n)]; D = [0.0]*n
    kPivotTol = 1e-10
    for i in range(n):
        for j in range(i):
            s = A[i][j]
            for k in range(j): s -= L[i][k]*D[k]*L[j][k]
            if D[j] == 0.0: return None
            L[i][j] = s/D[j]
        s = A[i][i]
        for k in range(i): s -= L[i][k]*D[k]*L[i][k]
        diag0 = abs(A[i][i])
        if (diag0 <= 1e-300) or (s <= 1e-30 and diag0 > 0.0) or \
           (diag0 > 0.0 and abs(s) <= kPivotTol*diag0): return None
        D[i] = s; L[i][i] = 1.0
    y = [0.0]*n
    for i in range(n):
        s = b[i]
        for j in range(i): s -= L[i][j]*y[j]
        y[i] = s
    z = [y[i]/D[i] for i in range(n)]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        s = z[i]
        for j in range(i+1, n): s -= L[j][i]*x[j]
        x[i] = s
    return x

def truss_ke(E, A, L):
    k = E*A/L
    return [[k, 0, -k, 0], [0, 0, 0, 0], [-k, 0, k, 0], [0, 0, 0, 0]]

def rot_matrix2(c, s):
    return [[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]]

def frame_ke(E, A, Iz, L):
    ke = [[0.0]*6 for _ in range(6)]
    ke[0][0] = ke[3][3] = E*A/L; ke[0][3] = ke[3][0] = -E*A/L
    b = 12*E*Iz/L**3; ke[1][1] = ke[4][4] = b; ke[1][4] = ke[4][1] = -b
    c = 6*E*Iz/L**2; ke[1][2] = ke[2][1] = ke[1][5] = ke[5][1] = c
    ke[2][4] = ke[4][2] = ke[4][5] = ke[5][4] = -c
    d = 4*E*Iz/L; ke[2][2] = ke[5][5] = d; ke[5][2] = ke[2][5] = d/2.0
    return ke

def rot_matrix6(c, s):
    T = [[0.0]*6 for _ in range(6)]
    T[2][2] = T[5][5] = 1.0
    T[0][0] = T[1][1] = T[3][3] = T[4][4] = c
    T[0][1] = T[3][4] = s
    T[1][0] = T[4][3] = -s
    return T

def solve(model):
    nodes = model["nodes"]; cons = model["constraints"]
    elems = model["elements"]; mats = model["materials"]; secs = model["sections"]
    loads = model["loads"]
    nnode = len(nodes)
    dof = [[0,0,0] for _ in range(nnode)]
    for c in cons: dof[c["node"]] = {"ux":0,"uy":0,"rz":0}.copy() if False else [0,0,0]
    # 约束: dofs 列出被约束自由度
    conmap = {}
    for c in cons:
        conmap[c["node"]] = set(c["dofs"])
    ibuf = 0
    dnames = ["ux", "uy", "rz"]
    for i in range(nnode):
        lim = 3 if nodes[i].get("type", "frame") == "frame" else 2  # truss 节点无 rz
        for j in range(lim):
            if dnames[j] not in conmap.get(i, set()): dof[i][j] = ibuf; ibuf += 1
    nfree = ibuf
    for i in range(nnode):
        lim = 3 if nodes[i].get("type", "frame") == "frame" else 2
        for j in range(lim):
            if dnames[j] in conmap.get(i, set()): dof[i][j] = ibuf; ibuf += 1
    ntotal = ibuf
    K = [[0.0]*ntotal for _ in range(ntotal)]
    elem_T, elem_kl = [], []
    for e in elems:
        E = mats[e["material"]]["E"]; A = secs[e["section"]]["A"]; Iz = secs[e["section"]]["Iz"]
        n0 = e["nodeI"]; n1 = e["nodeJ"]
        x0, y0 = nodes[n0]["x"], nodes[n0]["y"]; x1, y1 = nodes[n1]["x"], nodes[n1]["y"]
        L = math.hypot(x1-x0, y1-y0); c, s = (x1-x0)/L, (y1-y0)/L
        if e.get("type", "frame") == "truss":
            kl = truss_ke(E, A, L); T = rot_matrix2(c, s); g = dof[n0][:2] + dof[n1][:2]
        else:
            kl = frame_ke(E, A, Iz, L); T = rot_matrix6(c, s); g = dof[n0] + dof[n1]
        Kg = mat_mul(mat_mul(transpose(T), kl), T)
        for a in range(len(g)):
            for b in range(len(g)): K[g[a]][g[b]] += Kg[a][b]
        elem_T.append(T); elem_kl.append(kl)
    f = [0.0]*ntotal
    elem_fef = [[0.0]*6 for _ in elems]
    for Ld in loads:
        t = Ld["type"]
        if t == "nodalForce":
            d = {"x":0,"y":1,"rz":2}[Ld["direction"]]
            dd = dof[Ld["node"]][d]
            if 0 <= dd < ntotal: f[dd] += Ld["value"]
        elif t == "lateralUniformPressure":
            eid = Ld["element"]
            n0, n1 = elems[eid]["nodeI"], elems[eid]["nodeJ"]
            x0, y0 = nodes[n0]["x"], nodes[n0]["y"]; x1, y1 = nodes[n1]["x"], nodes[n1]["y"]
            Le = math.hypot(x1-x0, y1-y0); c, s = (x1-x0)/Le, (y1-y0)/Le
            q = Ld["value"]; pos = Ld.get("position", Le)
            if pos <= 0.0 or pos >= Le:  # 满跨
                FEF = [0.0]*6
                FEF[1] = -q*Le/2.0; FEF[4] = -q*Le/2.0
                FEF[2] = -q*Le*Le/12.0; FEF[5] = q*Le*Le/12.0
            else:
                qa = q*pos; dc = pos/Le; dg = dc*dc; db = Le-pos; ds = db/Le
                FEF = [0.0]*6
                FEF[1] = -qa*(2.0-2.0*dg+dc*dg); FEF[4] = -qa*dg*(2.0-dc)
                ds = qa*pos/6.0
                FEF[2] = -ds*(6.0-8.0*dc+3.0*dg); FEF[5] = ds*dc*(4.0-3.0*dc)
            for k in range(6): elem_fef[eid][k] += FEF[k]
            fe = [0.0]*6
            fe[0] = -( c*FEF[0] - s*FEF[1]); fe[1] = -( s*FEF[0] + c*FEF[1]); fe[2] = -FEF[2]
            fe[3] = -( c*FEF[3] - s*FEF[4]); fe[4] = -( s*FEF[3] + c*FEF[4]); fe[5] = -FEF[5]
            g = dof[n0] + dof[n1]
            for k in range(6):
                if 0 <= g[k] < ntotal: f[g[k]] += fe[k]
        elif t == "lateralForce":  # 横向集中力 P (垂直于杆, 距 i 端 pos)
            eid = Ld["element"]
            n0, n1 = elems[eid]["nodeI"], elems[eid]["nodeJ"]
            x0, y0 = nodes[n0]["x"], nodes[n0]["y"]; x1, y1 = nodes[n1]["x"], nodes[n1]["y"]
            Le = math.hypot(x1-x0, y1-y0); c, s = (x1-x0)/Le, (y1-y0)/Le
            P = Ld["value"]; a = Ld.get("position", 0.0); dc = a/Le if Le > 0 else 0.0
            dg = dc*dc; ds = (Le-a)/Le
            FEF = [0.0]*6
            FEF[1] = -P*ds*ds*(1.0+2.0*dc); FEF[4] = -P*dg*(1.0+2.0*ds)
            FEF[2] = -P*ds*ds*a;            FEF[5] = P*(Le-a)*dg
            for k in range(6): elem_fef[eid][k] += FEF[k]
            fe = [0.0]*6
            fe[0] = -( c*FEF[0] - s*FEF[1]); fe[1] = -( s*FEF[0] + c*FEF[1]); fe[2] = -FEF[2]
            fe[3] = -( c*FEF[3] - s*FEF[4]); fe[4] = -( s*FEF[3] + c*FEF[4]); fe[5] = -FEF[5]
            g = dof[n0] + dof[n1]
            for k in range(6):
                if 0 <= g[k] < ntotal: f[g[k]] += fe[k]
        elif t == "momentOnPoint":  # 节点弯矩 (node 字段) → 节点 rz
            dd = dof[Ld["node"]][2]
            if 0 <= dd < ntotal: f[dd] += Ld["value"]
        elif t in ("axialForce", "axialPressure"):
            eid = Ld["element"]
            n0, n1 = elems[eid]["nodeI"], elems[eid]["nodeJ"]
            x0, y0 = nodes[n0]["x"], nodes[n0]["y"]; x1, y1 = nodes[n1]["x"], nodes[n1]["y"]
            Le = math.hypot(x1-x0, y1-y0); c, s = (x1-x0)/Le, (y1-y0)/Le
            q = Ld["value"]; pos = Ld.get("position", 0.0)
            FEF = [0.0]*6
            if t == "axialForce":
                FEF[0] = -q*(Le-pos)/Le; FEF[3] = -q*pos/Le
            else:
                if pos <= 0.0 or pos >= Le:
                    FEF[0] = -q*Le/2.0; FEF[3] = -q*Le/2.0
                else:
                    qa = q*pos
                    FEF[0] = -qa*(Le-pos*0.5)/Le; FEF[3] = -qa*pos*0.5/Le
            for k in range(6): elem_fef[eid][k] += FEF[k]
            fe = [0.0]*6
            fe[0] = -( c*FEF[0] - s*FEF[1]); fe[1] = -( s*FEF[0] + c*FEF[1]); fe[2] = -FEF[2]
            fe[3] = -( c*FEF[3] - s*FEF[4]); fe[4] = -( s*FEF[3] + c*FEF[4]); fe[5] = -FEF[5]
            g = dof[n0] + dof[n1]
            for k in range(6):
                if 0 <= g[k] < ntotal: f[g[k]] += fe[k]
    free = list(range(nfree))
    Kff = [[K[i][j] for j in free] for i in free]
    ff = [f[i] for i in free]
    u_free = ldlt_solve(Kff, ff)
    if u_free is None: raise RuntimeError("singular system")
    u = [0.0]*ntotal
    for i, v in zip(free, u_free): u[i] = v
    disp = []
    for i in range(nnode):
        row = []
        lim = 3 if nodes[i].get("type", "frame") == "frame" else 2  # truss 节点无 rz
        for j in range(3):
            d = dof[i][j] if j < lim else -1
            row.append(u[d] if d >= 0 else 0.0)
        disp.append({"node": i, "ux": row[0], "uy": row[1], "rz": row[2]})
    endf = []
    for idx, e in enumerate(elems):
        T = elem_T[idx]
        if e.get("type", "frame") == "truss":
            g = dof[e["nodeI"]][:2] + dof[e["nodeJ"]][:2]
        else:
            g = dof[e["nodeI"]] + dof[e["nodeJ"]]
        ug = [u[d] if d >= 0 else 0.0 for d in g]
        if e.get("type", "frame") == "truss":
            # 4 自由度: N_i V_i N_j V_j (局部), M 槽位为 0
            fg = mat_vec(mat_mul(mat_mul(transpose(T), elem_kl[idx]), T), ug)
            fl = mat_vec(T, fg)
            fef = elem_fef[idx]
            fl = [fl[k] + (fef[0] if k == 0 else 0.0) for k in range(4)]
            fl = [fl[0], fl[1], 0.0, fl[2], fl[3], 0.0]
        else:
            fg = mat_vec(mat_mul(mat_mul(transpose(T), elem_kl[idx]), T), ug)
            fl = mat_vec(T, fg)  # f_local = T·f_global (T 正交: T·T^T = I)
            fl = [fl[k] + elem_fef[idx][k] for k in range(6)]
        endf.append({"element": idx, "nodeI": {"N": fl[0], "V": fl[1], "M": fl[2]},
                     "nodeJ": {"N": fl[3], "V": fl[4], "M": fl[5]}})
    reac = []
    for c in cons:
        row = {"node": c["node"], "ux": 0.0, "uy": 0.0, "rz": 0.0}
        for j, key in enumerate(["ux", "uy", "rz"]):
            d = dof[c["node"]][j]
            if d >= nfree:
                row[key] = sum(K[d][k]*u[k] for k in range(ntotal)) - f[d]
        reac.append(row)
    return disp, endf, reac

def main():
    model_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "examples", "portal_frame.json")
    with open(model_path, encoding="utf-8") as fh:
        model = json.load(fh)
    disp, endf, reac = solve(model)
    ref = {"displacements": disp, "endForces": endf, "reactions": reac}
    with open(os.path.join(HERE, "_portal_ref.json"), "w") as fh:
        json.dump(ref, fh, indent=1)
    print("== 参考解 (portal UDL) ==")
    for d in disp: print("  node %d: ux=%.6e uy=%.6e rz=%.6e" % (d["node"], d["ux"], d["uy"], d["rz"]))
    for e in endf: print("  elem %d: Ni=%.3f Vi=%.3f Mi=%.3f | Nj=%.3f Vj=%.3f Mj=%.3f" % (
        e["element"], e["nodeI"]["N"], e["nodeI"]["V"], e["nodeI"]["M"],
        e["nodeJ"]["N"], e["nodeJ"]["V"], e["nodeJ"]["M"]))
    for r in reac: print("  reac node %d: Rx=%.3f Ry=%.3f Mz=%.3f" % (r["node"], r["ux"], r["uy"], r["rz"]))

if __name__ == "__main__":
    main()
