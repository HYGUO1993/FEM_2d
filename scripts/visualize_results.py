import argparse
import re
import math
import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

def parse_results(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f.readlines()]
    idx = 0
    n_nodes = 0
    nodes = []
    elements = []
    displ = {}
    end_forces = {}
    reactions = {}
    while idx < len(lines):
        line = lines[idx].strip()
        if line.startswith("Total nodes:"):
            m = re.search(r"Total nodes:\s+(\d+)", line)
            if m:
                n_nodes = int(m.group(1))
        if line.startswith("Node data:"):
            idx += 2
            count = 0
            while idx < len(lines) and count < n_nodes:
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    continue
                parts = t.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    node_type = int(parts[0])
                    x = float(parts[1])
                    y = float(parts[2])
                    nodes.append((node_type, x, y))
                    count += 1
                idx += 1
            continue
        if line.startswith("Element data:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                parts = t.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    etype = int(parts[0])
                    i0 = int(parts[1])
                    i1 = int(parts[2])
                    sec = int(parts[3])
                    mat = int(parts[4])
                    elements.append((etype, i0, i1, sec, mat))
                idx += 1
            continue
        if line.startswith("Node Displacements:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                nums = re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', t, flags=re.IGNORECASE)
                if len(nums) >= 4:
                    i = int(float(nums[0]))
                    ux = float(nums[1])
                    uy = float(nums[2])
                    rz = float(nums[3])
                    displ[i] = (ux, uy, rz)
                idx += 1
            continue
        if line.startswith("Element End Forces:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                nums = re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', t, flags=re.IGNORECASE)
                if len(nums) >= 7:
                    e = int(float(nums[0]))
                    vals = list(map(float, nums[1:7]))
                    end_forces[e] = vals
                idx += 1
            continue
        if line.startswith("Support Reactions:"):
            idx += 2
            while idx < len(lines):
                t = lines[idx].strip()
                if not t:
                    idx += 1
                    break
                nums = re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', t, flags=re.IGNORECASE)
                if len(nums) >= 4:
                    node = int(float(nums[0]))
                    rx = float(nums[1])
                    ry = float(nums[2])
                    rz = float(nums[3])
                    reactions[node] = (rx, ry, rz)
                idx += 1
            continue
        idx += 1
    return nodes, elements, displ, end_forces, reactions

def plot_structure(nodes, elements, displ, end_forces, reactions, scale, out_path, show_reactions=True, reac_scale=None, show_element_force=True):
    xs0 = [n[1] for n in nodes]
    ys0 = [n[2] for n in nodes]
    spanx = (max(xs0) - min(xs0)) if xs0 else 1.0
    spany = (max(ys0) - min(ys0)) if ys0 else 1.0
    base = max(spanx, spany)
    if base <= 0:
        base = 1.0
    cap = 1000.0 * base
    xs1 = []
    ys1 = []
    for i in range(len(nodes)):
        ux, uy, rz = displ.get(i, (0.0, 0.0, 0.0))
        dx = scale * ux
        dy = scale * uy
        if not math.isfinite(dx) or abs(dx) > cap:
            dx = 0.0
        if not math.isfinite(dy) or abs(dy) > cap:
            dy = 0.0
        xs1.append(nodes[i][1] + dx)
        ys1.append(nodes[i][2] + dy)
    fig, ax = plt.subplots(figsize=(7, 7), dpi=130)
    # element force colormap
    axial_vals = []
    if show_element_force and end_forces:
        for idx_e, (etype, i0, i1, sec, mat) in enumerate(elements):
            ef = end_forces.get(idx_e)
            if ef is None:
                axial_vals.append(0.0)
                continue
            x0, y0 = xs0[i0], ys0[i0]
            x1, y1 = xs0[i1], ys0[i1]
            dx0 = x1 - x0
            dy0 = y1 - y0
            L = math.hypot(dx0, dy0)
            if L <= 0:
                axial_vals.append(0.0)
                continue
            ux = dx0 / L
            uy = dy0 / L
            fx_i, fy_i = ef[0], ef[1]
            fx_j, fy_j = ef[3], ef[4]
            n_i = fx_i * ux + fy_i * uy
            n_j = fx_j * ux + fy_j * uy
            axial = 0.5 * (n_i + n_j)
            axial_vals.append(abs(axial))
    vmin = min(axial_vals) if axial_vals else 0.0
    vmax = max(axial_vals) if axial_vals else 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("plasma")
    for etype, i0, i1, sec, mat in elements:
        ax.plot([xs0[i0], xs0[i1]], [ys0[i0], ys0[i1]], color="#999999", linewidth=2, label="_nolegend_")
    ax.scatter(xs0, ys0, s=30, color="#666666", label="original")
    for idx_e, (etype, i0, i1, sec, mat) in enumerate(elements):
        color = "#d62728"
        if show_element_force and axial_vals:
            color = cmap(norm(axial_vals[idx_e]))
        ax.plot([xs1[i0], xs1[i1]], [ys1[i0], ys1[i1]], color=color, linewidth=2, label="_nolegend_")
    ax.scatter(xs1, ys1, s=30, color="#d62728", label="deformed")
    if show_reactions and reactions:
        if reac_scale is None:
            reac_scale = 0.05 * base
        for node, (rx, ry, rz) in reactions.items():
            x = xs0[node]
            y = ys0[node]
            ax.quiver([x], [y], [rx * reac_scale], [ry * reac_scale], angles='xy', scale_units='xy', scale=1.0, color="#1f77b4", width=0.004, label="_nolegend_")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="upper right")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ttl = "Structure and Deformed Shape (scale={})".format(scale)
    if show_element_force and axial_vals:
        ttl += " | element color ~ axial"
    ax.set_title(ttl)
    if out_path:
        d = os.path.dirname(out_path)
        if d and not os.path.exists(d):
            os.makedirs(d)
        plt.savefig(out_path, bbox_inches="tight")
    else:
        plt.savefig("plot.png", bbox_inches="tight")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default="build/Results.dat")
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--out", default="build/plot.png")
    p.add_argument("--hide-reactions", action="store_true")
    p.add_argument("--hide-forces", action="store_true")
    p.add_argument("--reac-scale", type=float, default=None)
    args = p.parse_args()
    nodes, elements, displ, end_forces, reactions = parse_results(args.results)
    if not nodes or not elements:
        print("failed to parse results", file=sys.stderr)
        sys.exit(1)
    plot_structure(nodes, elements, displ, end_forces, reactions, args.scale, args.out, show_reactions=not args.hide_reactions, reac_scale=args.reac_scale, show_element_force=not args.hide_forces)
    print("saved", args.out)

if __name__ == "__main__":
    main()
