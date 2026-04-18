import math
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import colormaps
import matplotlib.patches as patches

from .fem_parser import FEMModelData
from .theme import PALETTE

class FEMVisualizer:
    """Handles matplotlib rendering for the FEM model."""
    
    def __init__(self, figure, canvas):
        self.figure = figure
        self.canvas = canvas
        self.ax = self.figure.add_subplot(111)
        
    def draw_empty(self):
        """Draw placeholder visualization."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title("Visualization Preview", color=PALETTE['text'])
        self.ax.text(
            0.5, 0.5,
            "1. Open a model file (.txt)\\n2. Click Solve (F5)\\n3. Explore the plot with pan/zoom tools",
            ha="center", va="center",
            transform=self.ax.transAxes,
            color=PALETTE['subtext0'],
            fontsize=11
        )
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        # Set colors directly on axes instead of relying on rcParams sometimes failing locally
        self.ax.set_facecolor(PALETTE['base'])
        self.figure.patch.set_facecolor(PALETTE['crust'])
        self.canvas.draw_idle()

    def _draw_constraint(self, x, y, cx, cy, cr, scale_base):
        """Draws constraint symbols (triangles, lines, etc.) based on fixed DOFs."""
        # scale_base is used to size the symbols appropriately
        s = scale_base * 0.05
        
        # fully fixed
        if cx < 0 and cy < 0 and cr < 0:
            # Draw a block hash (Horizontal base)
            self.ax.plot([x-s, x+s], [y, y], color=PALETTE['sapphire'], linewidth=3, zorder=4)
            for i in range(-2, 3):
                self.ax.plot([x+i*s/2.5, x+i*s/2.5-s/2], [y, y-s/2], color=PALETTE['sapphire'], linewidth=1.5, zorder=4)
                
        # pinned (cx, cy constrained, cr free)
        elif cx < 0 and cy < 0:
            triangle = patches.Polygon([[x, y], [x-s/1.5, y-s], [x+s/1.5, y-s]], 
                                       closed=True, facecolor=PALETTE['crust'], 
                                       edgecolor=PALETTE['sapphire'], zorder=4, linewidth=1.5)
            self.ax.add_patch(triangle)
            self.ax.plot([x-s, x+s], [y-s, y-s], color=PALETTE['sapphire'], linewidth=2, zorder=4)
            
        # roller (cy constrained only)
        elif cy < 0 and cx >= 0:
            triangle = patches.Polygon([[x, y], [x-s/1.5, y-s*0.8], [x+s/1.5, y-s*0.8]], 
                                       closed=True, facecolor=PALETTE['crust'], 
                                       edgecolor=PALETTE['sapphire'], zorder=4, linewidth=1.5)
            self.ax.add_patch(triangle)
            self.ax.plot([x-s, x+s], [y-s*1.1, y-s*1.1], color=PALETTE['sapphire'], linewidth=2, zorder=4)
            circ = patches.Circle((x-s/2, y-s*0.95), s*0.15, facecolor='none', edgecolor=PALETTE['sapphire'], zorder=4)
            circ2 = patches.Circle((x+s/2, y-s*0.95), s*0.15, facecolor='none', edgecolor=PALETTE['sapphire'], zorder=4)
            self.ax.add_patch(circ)
            self.ax.add_patch(circ2)
            
        # other partial constraints
        else:
            if cx < 0:
                self.ax.plot([x-s, x-s], [y-s, y+s], color=PALETTE['sapphire'], linewidth=2, zorder=4)
            if cy < 0:
                self.ax.plot([x-s, x+s], [y-s, y-s], color=PALETTE['sapphire'], linewidth=2, zorder=4)

    def render(self, model: FEMModelData, scale: float, reaction_scale: float,
               show_deformed: bool, show_forces: bool, show_reactions: bool, show_loads: bool = True,
               show_nodes: bool = False, show_elems: bool = False, show_bmd: bool = False):
        """Render the full model based on provided options."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.set_facecolor(PALETTE['base'])
        self.figure.patch.set_facecolor(PALETTE['crust'])

        if not model.nodes or not model.elements:
            self.draw_empty()
            return

        xs0 = [n[1] for n in model.nodes]
        ys0 = [n[2] for n in model.nodes]
        spanx = (max(xs0) - min(xs0)) if xs0 else 1.0
        spany = (max(ys0) - min(ys0)) if ys0 else 1.0
        base = max(spanx, spany)
        if base <= 0: base = 1.0
        cap = 1000.0 * base

        # Original Structure and Element Numbers
        for idx_e, (_, i0, i1, _, _) in enumerate(model.elements):
            x_m, y_m = (xs0[i0] + xs0[i1])/2, (ys0[i0] + ys0[i1])/2
            self.ax.plot([xs0[i0], xs0[i1]], [ys0[i0], ys0[i1]], color=PALETTE['overlay0'], linewidth=2.0, alpha=0.9, zorder=1)
            if show_elems:
                self.ax.text(x_m, y_m, f"E{idx_e}", color=PALETTE['blue'], fontsize=8, zorder=9,
                             ha='center', va='center',
                             bbox=dict(facecolor=PALETTE['crust'], alpha=0.6, edgecolor=PALETTE['surface1'], pad=1.0))
        
        # Nodes
        self.ax.scatter(xs0, ys0, s=30, color=PALETTE['overlay1'], label="Original", zorder=2)
        
        # Node Numbers
        if show_nodes:
            for i, (x, y) in enumerate(zip(xs0, ys0)):
                self.ax.text(x + base*0.015, y + base*0.015, str(i), color=PALETTE['text'], 
                             fontsize=8, zorder=10, alpha=0.9,
                             bbox=dict(facecolor=PALETTE['mantle'], alpha=0.7, edgecolor='none', pad=1))

        # Constraints
        for node_idx, (cx, cy, cr) in model.constraints.items():
            if node_idx < len(xs0):
                self._draw_constraint(xs0[node_idx], ys0[node_idx], cx, cy, cr, base)

        # External Loads (pre-solve visual indicator)
        if show_loads and not model.is_solved:
            # Draw external loads
            max_load = max([abs(l[2]) for l in model.loads] + [1e-6])
            for type_load, dir_load, val, elem, node, _, _, _ in model.loads:
                if type_load == 1:
                    # Node load
                    load_scale = (base * 0.15) / max_load
                    x, y = xs0[node], ys0[node]
                    if dir_load == 0:
                        self.ax.arrow(x - val * load_scale, y, val * load_scale, 0, head_width=base*0.02, head_length=base*0.03,
                                      fc=PALETTE['peach'], ec=PALETTE['peach'], zorder=4, length_includes_head=True)
                        self.ax.text(x - val * load_scale, y, f"{val:.1f}", color=PALETTE['peach'], fontsize=8)
                    elif dir_load == 1:
                        self.ax.arrow(x, y - val * load_scale, 0, val * load_scale, head_width=base*0.02, head_length=base*0.03,
                                      fc=PALETTE['peach'], ec=PALETTE['peach'], zorder=4, length_includes_head=True)
                        self.ax.text(x, y - val * load_scale, f"{val:.1f}", color=PALETTE['peach'], fontsize=8)

        # Deformed Shape & Results
        if model.is_solved and show_deformed:
            xs1, ys1 = [], []
            for i in range(len(model.nodes)):
                ux, uy, _ = model.displacements.get(i, (0.0, 0.0, 0.0))
                dx = scale * ux
                dy = scale * uy
                if not math.isfinite(dx) or abs(dx) > cap: dx = 0.0
                if not math.isfinite(dy) or abs(dy) > cap: dy = 0.0
                xs1.append(xs0[i] + dx)
                ys1.append(ys0[i] + dy)

            # Axial Forces calculation for color coding
            axial_vals = []
            for idx_e, (_, i0, i1, _, _) in enumerate(model.elements):
                ef = model.end_forces.get(idx_e)
                if ef is None:
                    axial_vals.append(0.0)
                    continue
                dx0, dy0 = xs0[i1] - xs0[i0], ys0[i1] - ys0[i0]
                length = math.hypot(dx0, dy0)
                if length <= 0:
                    axial_vals.append(0.0)
                    continue
                vec_x, vec_y = dx0/length, dy0/length
                n_i = ef[0] * vec_x + ef[1] * vec_y
                n_j = ef[3] * vec_x + ef[4] * vec_y
                axial_vals.append(abs(0.5*(n_i + n_j)))
                
            vmin = min(axial_vals) if axial_vals else 0.0
            vmax = max(axial_vals) if axial_vals else 1.0
            if vmax <= vmin: vmax = vmin + 1.0
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
            cmap = colormaps["plasma"]

            for idx_e, (_, i0, i1, _, _) in enumerate(model.elements):
                color = PALETTE['red']
                if show_forces and axial_vals:
                    color = cmap(norm(axial_vals[idx_e]))
                self.ax.plot([xs1[i0], xs1[i1]], [ys1[i0], ys1[i1]], color=color, linewidth=2.5, zorder=6)
                
            self.ax.scatter(xs1, ys1, s=35, color=PALETTE['red'], label="Deformed", zorder=7)

            # Reactions
            if show_reactions and model.reactions:
                max_r = max([abs(r[0]) for r in model.reactions.values()] + [abs(r[1]) for r in model.reactions.values()] + [1e-6])
                arrow_scale = (base * 0.15) / max_r
                for node, (rx, ry, _) in model.reactions.items():
                    if node < len(xs0):
                        x, y = xs0[node], ys0[node]
                        if abs(rx) > 1e-3:
                            self.ax.arrow(x - rx * arrow_scale, y, rx * arrow_scale, 0, head_width=base*0.02, head_length=base*0.03, length_includes_head=True,
                                          fc=PALETTE['maroon'], ec=PALETTE['maroon'], zorder=8)
                            self.ax.text(x - rx * arrow_scale, y, f"{rx:.1f}", color=PALETTE['maroon'], fontsize=7,
                                         bbox=dict(facecolor=PALETTE['base'], alpha=0.5, edgecolor='none', pad=0.5))
                        if abs(ry) > 1e-3:
                            self.ax.arrow(x, y - ry * arrow_scale, 0, ry * arrow_scale, head_width=base*0.02, head_length=base*0.03, length_includes_head=True,
                                          fc=PALETTE['maroon'], ec=PALETTE['maroon'], zorder=8)
                            self.ax.text(x, y - ry * arrow_scale, f"{ry:.1f}", color=PALETTE['maroon'], fontsize=7,
                                         bbox=dict(facecolor=PALETTE['base'], alpha=0.5, edgecolor='none', pad=0.5))
                                         
            # BMD (Bending Moment Diagram)
            if show_bmd and model.end_forces:
                max_m = 0.0
                for idx_e, (t, i0, i1, _, _) in enumerate(model.elements):
                    if t == 2 and idx_e in model.end_forces:
                        ef = model.end_forces[idx_e]
                        max_m = max(max_m, abs(ef[2]), abs(ef[5]))
                
                if max_m > 1e-6:
                    bmd_scale = (base * 0.15) / max_m
                    for idx_e, (t, i0, i1, _, _) in enumerate(model.elements):
                        if t == 2 and idx_e in model.end_forces:
                            ef = model.end_forces[idx_e]
                            mzi, mzj = ef[2], ef[5]
                            # BMD should be drawn on the original structure (xs0, ys0) to avoid deformation distortion
                            x_i, y_i = xs0[i0], ys0[i0]
                            x_j, y_j = xs0[i1], ys0[i1]
                            dx, dy = x_j - x_i, y_j - y_i
                            L = math.hypot(dx, dy)
                            if L > 0:
                                nx, ny = -dy/L, dx/L
                                # Continuous moment convention
                                val_i = -mzi * bmd_scale
                                val_j = mzj * bmd_scale 
                                vx_i, vy_i = nx * val_i, ny * val_i
                                vx_j, vy_j = nx * val_j, ny * val_j 
                                poly_x = [x_i, x_i + vx_i, x_j + vx_j, x_j]
                                poly_y = [y_i, y_i + vy_i, y_j + vy_j, y_j]
                                self.ax.fill(poly_x, poly_y, color=PALETTE['green'], alpha=0.3, zorder=5)
                                self.ax.plot([x_i + vx_i, x_j + vx_j], [y_i + vy_i, y_j + vy_j], color=PALETTE['green'], linewidth=1.5, zorder=6)
                                # Labels for max moments
                                if abs(mzi) > 1e-2:
                                    self.ax.text(x_i + vx_i*1.1, y_i + vy_i*1.1, f"{-mzi:.1f}", color=PALETTE['green'], fontsize=7,
                                                 bbox=dict(facecolor=PALETTE['base'], alpha=0.5, edgecolor='none', pad=0.5))
                                if abs(mzj) > 1e-2:
                                    self.ax.text(x_j + vx_j*1.1, y_j + vy_j*1.1, f"{mzj:.1f}", color=PALETTE['green'], fontsize=7,
                                                 bbox=dict(facecolor=PALETTE['base'], alpha=0.5, edgecolor='none', pad=0.5))

            # Colorbar
            if show_forces and axial_vals:
                sm = colormaps["plasma"]
                try:
                    from matplotlib.cm import ScalarMappable
                    mappable = ScalarMappable(norm=norm, cmap=sm)
                    mappable.set_array([])
                    cbar = self.figure.colorbar(mappable, ax=self.ax, fraction=0.045, pad=0.03)
                    cbar.set_label("|Axial Force|", color=PALETTE['text'])
                    cbar.ax.yaxis.set_tick_params(color=PALETTE['subtext0'], labelcolor=PALETTE['text'])
                    cbar.outline.set_edgecolor(PALETTE['surface1'])
                except Exception as e:
                    print("Colorbar error:", e)

        # Plot Settings
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.grid(True, linestyle="--", alpha=0.3, color=PALETTE['surface0'])
        self.ax.set_xlabel("X", color=PALETTE['text'], fontweight='bold')
        self.ax.set_ylabel("Y", color=PALETTE['text'], fontweight='bold')
        self.ax.tick_params(colors=PALETTE['subtext0'])
        for spine in self.ax.spines.values():
            spine.set_color(PALETTE['surface2'])
            
        title = f"FEM Analysis Visualization"
        if model.is_solved:
            title += f" (Disp. Scale: {scale:g})"
        self.ax.set_title(title, color=PALETTE['text'], pad=10, fontweight='bold')
        
        # Legend styling
        legend = self.ax.legend(loc="upper right", facecolor=PALETTE['mantle'], edgecolor=PALETTE['surface1'])
        for text in legend.get_texts():
            text.set_color(PALETTE['text'])

        self.canvas.draw_idle()
