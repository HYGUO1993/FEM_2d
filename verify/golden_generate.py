#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Golden value generator for FEM_2d regression tests.
Fully self-contained reference FEM solver (dense assembly + correct LDLT),
independent from the C++ implementation. Produces tests/golden/golden.json
with node displacements, element end forces (LOCAL coordinates: N/V/M),
and support reactions for the 3 canonical cases.

Usage: python golden_generate.py
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "tests", "golden", "golden.json")

# ---------------------------------------------------------------- helpers
def parse_input(path):
    with open(path) as f:
        toks = f.read().split()
    nNode, nCons, nElem, nMat, nSec, nLoad = (int(toks[i]) for i in range(6))
    idx = 6
    nodes, cons, elems, mats, secs, loads = [], [], [], [], [], []
    for _ in range(nNode):
        nodes.append((int(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
    for _ in range(nCons):
        cons.append((int(toks[idx]), int(toks[idx+1]), int(toks[idx+2]), int(toks[idx+3]))); idx += 4
    for _ in range(nElem):
        elems.append((int(toks[idx]), int(toks[idx+1]), int(toks[idx+2]), int(toks[idx+3]), 0)); idx += 4
    for _ in range(nMat):
        mats.append((float(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
    for _ in range(nSec):
        secs.append((float(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
    for _ in range(nLoad):
        loads.append(tuple(float(toks[idx+k]) for k in range(8))); idx += 8
    return nodes, cons, elems, mats, secs, loads

def dof_numbering(nodes, cons):
    n = len(nodes); dof = [[0, 0, 0] for _ in range(n)]
    for (cn, a, b, c) in cons:
        dof[cn] = [a, b, c]
    ibuf = 0
    for i in range(n):
        lim = 3 if nodes[i][0] == 2 else 2
        for j in range(lim):
            if dof[i][j] == 0: dof[i][j] = ibuf; ibuf += 1
    nfree = ibuf
    for i in range(n):
        lim = 3 if nodes[i][0] == 2 else 2
        for j in range(lim):
            if dof[i][j] == -1: dof[i][j] = ibuf; ibuf += 1
    return dof, nfree, ibuf

def rot_matrix2(c, s):
    """4x4 rotation for truss: x' = c*x + s*y, y' = -s*x + c*y (matches C++ TrussElemStiffCalcu)."""
    # block diag [[c, s],[-s, c]] twice
    return [[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]]

def rot_matrix6(c, s):
    """6x6 rotation for frame (x' axis along element)."""
    T = [[0.0]*6 for _ in range(6)]
    T[2][2] = T[5][5] = 1.0
    T[0][0] = T[1][1] = T[3][3] = T[4][4] = c
    T[0][1] = T[3][4] = s
    T[1][0] = T[4][3] = -s
    return T

def truss_ke(E, A, L):
    k = E*A/L
    return [[k, 0, -k, 0], [0, 0, 0, 0], [-k, 0, k, 0], [0, 0, 0, 0]]

def frame_ke(E, A, Iz, L):
    ke = [[0.0]*6 for _ in range(6)]
    ke[0][0] = ke[3][3] = E*A/L; ke[0][3] = ke[3][0] = -E*A/L
    b = 12*E*Iz/L**3; ke[1][1] = ke[4][4] = b; ke[1][4] = ke[4][1] = -b
    c = 6*E*Iz/L**2; ke[1][2] = ke[2][1] = ke[1][5] = ke[5][1] = c
    ke[2][4] = ke[4][2] = ke[4][5] = ke[5][4] = -c
    d = 4*E*Iz/L; ke[2][2] = ke[5][5] = d; ke[5][2] = ke[2][5] = d/2.0
    return ke

def mat_vec(A, x):
    return [sum(A[i][j]*x[j] for j in range(len(x))) for i in range(len(A))]

def mat_mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def transpose(A):
    return [list(r) for r in zip(*A)]

def ldlt_solve(A, b):
    """LDLT solve with singularity detection (mirrors barsystem.cpp LDLTSolve).
    Returns the solution vector, or None if the system is singular/ill-conditioned.
    """
    n = len(A); L = [[0.0]*n for _ in range(n)]; D = [0.0]*n
    kPivotTol = 1e-10
    for i in range(n):
        for j in range(i):
            s = A[i][j]
            for k in range(j): s -= L[i][k]*D[k]*L[j][k]
            if D[j] == 0.0:
                return None
            L[i][j] = s/D[j]
        s = A[i][i]
        for k in range(i): s -= L[i][k]*D[k]*L[i][k]
        diag0 = abs(A[i][i])
        if (diag0 <= 1e-300) or (s <= 1e-30 and diag0 > 0.0) or \
           (diag0 > 0.0 and abs(s) <= kPivotTol*diag0):
            return None
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

def solve_case(case_path):
    nodes, cons, elems, mats, secs, loads = parse_input(case_path)
    dof, nfree, ntotal = dof_numbering(nodes, cons)
    n = ntotal
    K = [[0.0]*n for _ in range(n)]
    # store per-element T and K_local for force recovery
    elem_T, elem_kl, elem_g = [], [], []
    for (itype, n0, n1, sec, mat) in elems:
        E, A, Iz = mats[mat][0], secs[sec][0], secs[sec][1]
        x0, y0 = nodes[n0][1], nodes[n0][2]; x1, y1 = nodes[n1][1], nodes[n1][2]
        L = math.hypot(x1-x0, y1-y0)
        c, s = (x1-x0)/L, (y1-y0)/L
        if itype == 1:
            kl = truss_ke(E, A, L); T = rot_matrix2(c, s); g = dof[n0][:2] + dof[n1][:2]
        else:
            kl = frame_ke(E, A, Iz, L); T = rot_matrix6(c, s); g = dof[n0] + dof[n1]
        Kg = mat_mul(mat_mul(transpose(T), kl), T)
        for a in range(len(g)):
            for b in range(len(g)):
                K[g[a]][g[b]] += Kg[a][b]
        elem_T.append(T); elem_kl.append(kl); elem_g.append(g)
    f = [0.0]*n
    # per-element fixed-end force accumulation (local), for end-force recovery
    elem_fef = [[0.0]*6 for _ in elems]
    for Ld in loads:
        typ, direct, val = int(Ld[0]), int(Ld[1]), Ld[2]
        ln = int(Ld[4])
        if typ == 1:  # FORCE_ON_NODE
            d = dof[ln][direct]
            if 0 <= d < n: f[d] += val
        elif typ == 3:  # LATERAL_UNIFORM_PRESSURE (full-span udl on element)
            eid = int(Ld[3])
            q = val
            a0, b0 = nodes[elems[eid][1]][1], nodes[elems[eid][1]][2]
            a1, b1 = nodes[elems[eid][2]][1], nodes[elems[eid][2]][2]
            Le = math.hypot(a1-a0, b1-b0)
            c, s = (a1-a0)/Le, (b1-b0)/Le
            # fixed-end forces (local, same sign convention as C++):
            FEF = [0.0]*6
            FEF[1] = -q*Le/2.0
            FEF[4] = -q*Le/2.0
            FEF[2] = -q*Le*Le/12.0
            FEF[5] =  q*Le*Le/12.0
            for k in range(6): elem_fef[eid][k] += FEF[k]
            # equivalent nodal load = -FEF, to global
            fe = [0.0]*6
            fe[0] = -( c*FEF[0] - s*FEF[1]); fe[1] = -( s*FEF[0] + c*FEF[1]); fe[2] = -FEF[2]
            fe[3] = -( c*FEF[3] - s*FEF[4]); fe[4] = -( s*FEF[3] + c*FEF[4]); fe[5] = -FEF[5]
            g = dof[elems[eid][1]] + dof[elems[eid][2]]
            for k in range(6):
                if 0 <= g[k] < n: f[g[k]] += fe[k]
    free = list(range(nfree))
    Kff = [[K[i][j] for j in free] for i in free]
    ff = [f[i] for i in free]
    u_free = ldlt_solve(Kff, ff)
    if u_free is None:
        raise RuntimeError("singular/ill-conditioned system for %s (mechanism or missing constraints)" % case_path)
    u = [0.0]*n
    for i, v in zip(free, u_free): u[i] = v
    # node displacements (3 dof per node, report 0 for constrained/absent)
    disp = []
    for i in range(len(nodes)):
        row = []
        for j in range(3):
            d = dof[i][j]
            row.append(u[d] if d >= 0 else 0.0)
        disp.append(row)
    # element end forces in LOCAL coordinates (N V M at each end)
    endf = []
    for idx, (itype, n0, n1, sec, mat) in enumerate(elems):
        T = elem_T[idx]; g = elem_g[idx]
        ug = [u[d] if d >= 0 else 0.0 for d in g]
        fg = mat_vec(mat_mul(mat_mul(transpose(T), elem_kl[idx]), T), ug)  # global end forces
        fl = mat_vec(transpose(T), fg)                                     # local end forces
        # add fixed-end forces: true end force = k*u + FEF
        fef = elem_fef[idx]
        if itype == 1:
            # truss local vector is 4-long (N_i V_i N_j V_j); FEF only 6-layout for frame
            n_fef = [fef[0], fef[1], fef[3], fef[4]]
            fl = [fl[k] + n_fef[k] for k in range(4)]
            endf.append([fl[0], fl[1], 0.0, fl[2], fl[3], 0.0])            # N_i V_i _ N_j V_j _
        else:
            fl = [fl[k] + fef[k] for k in range(6)]
            endf.append(fl)
    # reactions: R = K*u - f at constrained dofs
    reac = []
    for (cn, a, b, c) in cons:
        row = []
        for j in range(3):
            d = dof[cn][j]
            if d >= nfree:
                Ku = sum(K[d][k]*u[k] for k in range(n))
                row.append(Ku - f[d])
            else:
                row.append(0.0)
        reac.append(row)
    return nodes, disp, endf, reac

# ---------------------------------------------------------------- main
cases = [
    ("cantilever_frame", "tests/golden/inputs/test05.txt"),
    ("simply_supported_beam", "tests/golden/inputs/test_beam.txt"),
    ("triangle_truss", "tests/golden/inputs/test_triangle.txt"),
    ("cantilever_udl", "tests/golden/inputs/test_cantilever_udl.txt"),
]
result = {"cases": {}}
for name, fname in cases:
    nodes, disp, endf, reac = solve_case(os.path.join(ROOT, fname))
    result["cases"][name] = {
        "input": fname,
        "n_nodes": len(nodes),
        "displacements": disp,
        "end_forces": endf,
        "reactions": reac,
    }
with open(OUT, "w") as f:
    json.dump(result, f, indent=1)
print("Wrote", OUT)
for name in result["cases"]:
    c = result["cases"][name]
    print(f"\n== {name} ==")
    print("  disp:", c["displacements"])
    print("  endf:", c["end_forces"])
    print("  reac:", c["reactions"])
