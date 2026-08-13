#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent verification harness for FEM_2d.
Replicates the C++ element-stiffness + DOF numbering + skyline assembly using
dense matrices (no skyline bugs), then solves with BOTH:
  (A) a CORRECT LDL^T back-substitution  (reference / analytic)
  (B) the program's buggy back-substitution (uses pB[j] instead of z[j])
and compares against fem_run.exe's actual output.
Pure stdlib, no numpy required.
"""
import subprocess, os, re, math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, "build", "bin", "fem_run.exe")
INPUT = os.path.join(ROOT, "test05.txt")
OUT = os.path.join(HERE, "Results_frame_verify.dat")

# ---- parse test05.txt (cantilever frame) ----
with open(INPUT) as f:
    toks = f.read().split()
# header: 2 1 1 1 1 1
nNode, nCons, nElem, nMat, nSec, nLoad = (int(toks[0]), int(toks[1]), int(toks[2]),
                                          int(toks[3]), int(toks[4]), int(toks[5]))
idx = 6
nodes = []
for _ in range(nNode):
    nodes.append((int(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
cons = []
for _ in range(nCons):
    cons.append((int(toks[idx]), int(toks[idx+1]), int(toks[idx+2]), int(toks[idx+3]))); idx += 4
elems = []
for _ in range(nElem):
    itype = int(toks[idx]); n0 = int(toks[idx+1]); n1 = int(toks[idx+2]); sec = int(toks[idx+3])
    mat = 0
    # material optional on same line after section
    rest = toks[idx+4:idx+5]
    elems.append((itype, n0, n1, sec, mat)); idx += 4
mats = []
for _ in range(nMat):
    mats.append((float(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
secs = []
for _ in range(nSec):
    secs.append((float(toks[idx]), float(toks[idx+1]), float(toks[idx+2]))); idx += 3
loads = []
for _ in range(nLoad):
    loads.append(tuple(float(toks[idx+k]) for k in range(9))); idx += 9

# ---- DOF numbering (replicate DOFIndexCalcu) ----
TRUSS_NODE, FRAME_NODE = 1, 2
dof = [[0,0,0] for _ in range(nNode)]
for (cn, a, b, c) in cons:
    dof[cn] = [a, b, c]
ibuf = 0
for i in range(nNode):
    ntype = nodes[i][0]
    lim = 3 if ntype == FRAME_NODE else 2
    for j in range(lim):
        if dof[i][j] == 0:
            dof[i][j] = ibuf; ibuf += 1
nFree = ibuf
for i in range(nNode):
    ntype = nodes[i][0]
    lim = 3 if ntype == FRAME_NODE else 2
    for j in range(lim):
        if dof[i][j] == -1:
            dof[i][j] = ibuf; ibuf += 1
nTotal = ibuf
print("DOF numbering:", dof, "nFree=", nFree, "nTotal=", nTotal)

def frame_local_ke(E, A, Iz, L):
    """Standard 6x6 frame element stiffness in LOCAL coords (node i then j)."""
    ke = [[0.0]*6 for _ in range(6)]
    ea_l = E*A/L
    ei = E*Iz
    ke[0][0] = ke[3][3] = ea_l; ke[0][3] = ke[3][0] = -ea_l
    b = 12*ei/(L**3); ke[1][1] = ke[4][4] = b; ke[1][4] = ke[4][1] = -b
    c = 6*ei/(L**2); ke[1][2] = ke[2][1] = ke[1][5] = ke[5][1] = c
    ke[2][4] = ke[4][2] = ke[4][5] = ke[5][4] = -c
    d = 4*ei/L; ke[2][2] = ke[5][5] = d; ke[5][2] = ke[2][5] = d/2.0
    return ke

# assemble dense full 6x6 (global = local since element along x)
N = nTotal
K = [[0.0]*N for _ in range(N)]
for (itype, n0, n1, sec, mat) in elems:
    E, A, Iz = mats[mat][0], secs[sec][0], secs[sec][1]
    x0,y0 = nodes[n0][1], nodes[n0][2]; x1,y1 = nodes[n1][1], nodes[n1][2]
    L = math.hypot(x1-x0, y1-y0)
    ke = frame_local_ke(E, A, Iz, L)  # cos=1,sin=0 -> global==local
    g = dof[n0] + dof[n1]
    for a in range(6):
        for b in range(6):
            K[g[a]][g[b]] += ke[a][b]

# load vector (free DOFs)
f = [0.0]*N
for Ld in loads:
    typ, direct, val, _, ln, *_ = Ld
    if int(typ) == 1:  # FORCE_ON_NODE
        d = int(dof[int(ln)][int(direct)])
        if 0 <= d < N:
            f[d] += val

# reduced free-free system
free = list(range(nFree))
Kff = [[K[i][j] for j in free] for i in free]
ff = [f[i] for i in free]

def ldlt_correct(A, b):
    n = len(A); Lm = [[0.0]*n for _ in range(n)]; D = [0.0]*n
    for i in range(n):
        for j in range(i):
            s = A[i][j]
            for k in range(j): s -= Lm[i][k]*D[k]*Lm[j][k]
            Lm[i][j] = s/D[j] if D[j] != 0 else 0.0
        s = A[i][i]
        for k in range(i): s -= Lm[i][k]*D[k]*Lm[i][k]
        D[i] = s
        Lm[i][i] = 1.0
    y = [0.0]*n
    for i in range(n):
        s = b[i]
        for j in range(i): s -= Lm[i][j]*y[j]
        y[i] = s
    z = [y[i]/D[i] if D[i]!=0 else 0.0 for i in range(n)]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        s = z[i]
        for j in range(i+1, n): s -= Lm[j][i]*z[j]   # CORRECT: uses z[j]
        x[i] = s
    return x

def ldlt_buggy(A, b):
    """Replicates fem_run's LDLTSolve: back-sub uses pB[j] instead of z[j]."""
    n = len(A); Lm = [[0.0]*n for _ in range(n)]; D = [0.0]*n
    for i in range(n):
        for j in range(i):
            s = A[i][j]
            for k in range(j): s -= Lm[i][k]*D[k]*Lm[j][k]
            Lm[i][j] = s/D[j] if D[j] != 0 else 0.0
        s = A[i][i]
        for k in range(i): s -= Lm[i][k]*D[k]*Lm[i][k]
        D[i] = s
        Lm[i][i] = 1.0
    y = [0.0]*n
    for i in range(n):
        s = b[i]
        for j in range(i): s -= Lm[i][j]*y[j]
        y[i] = s
    z = [y[i]/D[i] if D[i]!=0 else 0.0 for i in range(n)]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        s = z[i]
        for j in range(i+1, n): s -= Lm[j][i]*b[j]   # BUG: uses b[j] (original RHS, now overwritten)
        x[i] = s
    return x

xc = ldlt_correct(Kff, ff)
xb = ldlt_buggy(Kff, ff)
print("\n--- Correct LDLT (reference) ---")
print("Uy(node1) =", xc[1], " Rz(node1) =", xc[2])
print("--- Buggy  LDLT (replicates fem_run) ---")
print("Uy(node1) =", xb[1], " Rz(node1) =", xb[2])

# analytic cantilever tip deflection
E, A, Iz = mats[0][0], secs[0][0], secs[0][1]
L = 2.0; P = 1000.0
defl = P*L**3/(3*E*Iz)
print("\nAnalytic cantilever: Uy =", -defl, " fixed-end M =", P*L)

# ---- run fem_run.exe and parse its results ----
r = subprocess.run([EXE, "--no-stiff", "--input", INPUT, "--output", OUT],
                    cwd=ROOT, capture_output=True, text=True)
with open(OUT) as fh:
    txt = fh.read()
m = re.search(r"Node\s+Ux\s+Uy\s+Rz\s+([\d\s\.\-]+)", txt)
# parse node displacements
disp_lines = []
for line in txt.splitlines():
    if re.match(r"^\s*\d+\s+[-0-9.]+", line) and "Node" not in line:
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit():
            disp_lines.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))
print("\nfem_run.exe Node Displacements:", disp_lines)
