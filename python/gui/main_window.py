"""
Main Window for FEM_2d GUI Application

This module implements the main application window with menu bar, toolbar,
status bar, and docked panels for model tree and properties.
"""

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QMenuBar, QMenu, QToolBar, QStatusBar, QDockWidget,
    QTreeWidget, QTreeWidgetItem, QTextEdit, QLabel,
    QFileDialog, QMessageBox, QSplitter, QPushButton,
    QCheckBox, QDoubleSpinBox, QFormLayout
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QAction, QIcon, QKeySequence

import os
import subprocess
import re
import math

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import colors as mcolors
from matplotlib import colormaps


class MainWindow(QMainWindow):
    """Main application window."""

    def __init__(self):
        super().__init__()

        # Window properties
        self.setWindowTitle("FEM_2d - 2D Finite Element Analysis")
        self.setGeometry(100, 100, 1400, 900)

        # Current model
        self.current_file = None
        self.model = None
        self.is_modified = False
        self.project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        self.last_results_path = os.path.join(self.project_root, "build", "Results_gui.dat")
        self.last_plot_path = os.path.join(self.project_root, "build", "plot_gui.png")
        self.results_data = None

        # Setup UI
        self._create_actions()
        self._create_menus()
        self._create_toolbars()
        self._create_status_bar()
        self._create_docks()
        self._create_central_widget()

        # Update UI state
        self._update_window_title()
        self._update_actions()

    def _create_actions(self):
        """Create application actions."""
        # File actions
        self.new_action = QAction("&New", self)
        self.new_action.setShortcut(QKeySequence.New)
        self.new_action.setStatusTip("Create a new model")
        self.new_action.triggered.connect(self.new_model)

        self.open_action = QAction("&Open...", self)
        self.open_action.setShortcut(QKeySequence.Open)
        self.open_action.setStatusTip("Open an existing model")
        self.open_action.triggered.connect(self.open_model)

        self.save_action = QAction("&Save", self)
        self.save_action.setShortcut(QKeySequence.Save)
        self.save_action.setStatusTip("Save the current model")
        self.save_action.triggered.connect(self.save_model)

        self.save_as_action = QAction("Save &As...", self)
        self.save_as_action.setShortcut(QKeySequence.SaveAs)
        self.save_as_action.setStatusTip("Save the model with a new name")
        self.save_as_action.triggered.connect(self.save_model_as)

        self.exit_action = QAction("E&xit", self)
        self.exit_action.setShortcut(QKeySequence.Quit)
        self.exit_action.setStatusTip("Exit the application")
        self.exit_action.triggered.connect(self.close)

        # Edit actions
        self.undo_action = QAction("&Undo", self)
        self.undo_action.setShortcut(QKeySequence.Undo)
        self.undo_action.setEnabled(False)

        self.redo_action = QAction("&Redo", self)
        self.redo_action.setShortcut(QKeySequence.Redo)
        self.redo_action.setEnabled(False)

        # Analysis actions
        self.solve_action = QAction("&Solve", self)
        self.solve_action.setShortcut(Qt.Key_F5)
        self.solve_action.setStatusTip("Run FEM analysis and generate visualization")
        self.solve_action.triggered.connect(self.solve_model)

        self.open_example_action = QAction("Open &Example (test05)", self)
        self.open_example_action.setStatusTip("Use test05.txt from project root")
        self.open_example_action.triggered.connect(self.open_default_example)

        # View actions
        self.show_deformed_action = QAction("Show &Deformed", self)
        self.show_deformed_action.setShortcut(Qt.Key_F6)
        self.show_deformed_action.setCheckable(True)
        self.show_deformed_action.setChecked(True)
        self.show_deformed_action.triggered.connect(self._refresh_visualization)

        self.show_forces_action = QAction("Show &Forces", self)
        self.show_forces_action.setShortcut(Qt.Key_F7)
        self.show_forces_action.setCheckable(True)
        self.show_forces_action.setChecked(True)
        self.show_forces_action.triggered.connect(self._refresh_visualization)

        self.show_reactions_action = QAction("Show &Reactions", self)
        self.show_reactions_action.setShortcut(Qt.Key_F8)
        self.show_reactions_action.setCheckable(True)
        self.show_reactions_action.setChecked(True)
        self.show_reactions_action.triggered.connect(self._refresh_visualization)

        self.export_plot_action = QAction("&Export Plot...", self)
        self.export_plot_action.setStatusTip("Export current visualization to image")
        self.export_plot_action.triggered.connect(self.export_plot)

        # Help actions
        self.about_action = QAction("&About", self)
        self.about_action.setStatusTip("About FEM_2d")
        self.about_action.triggered.connect(self.show_about)

    def _create_menus(self):
        """Create application menus."""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(self.new_action)
        file_menu.addAction(self.open_action)
        file_menu.addAction(self.open_example_action)
        file_menu.addSeparator()
        file_menu.addAction(self.save_action)
        file_menu.addAction(self.save_as_action)
        file_menu.addSeparator()
        file_menu.addAction(self.exit_action)

        # Edit menu
        edit_menu = menubar.addMenu("&Edit")
        edit_menu.addAction(self.undo_action)
        edit_menu.addAction(self.redo_action)

        # Analysis menu
        analysis_menu = menubar.addMenu("&Analysis")
        analysis_menu.addAction(self.solve_action)
        analysis_menu.addAction(self.export_plot_action)

        # View menu
        view_menu = menubar.addMenu("&View")
        view_menu.addAction(self.show_deformed_action)
        view_menu.addAction(self.show_forces_action)
        view_menu.addAction(self.show_reactions_action)

        # Help menu
        help_menu = menubar.addMenu("&Help")
        help_menu.addAction(self.about_action)

    def _create_toolbars(self):
        """Create application toolbars."""
        # Main toolbar
        main_toolbar = self.addToolBar("Main")
        main_toolbar.setMovable(False)
        main_toolbar.addAction(self.new_action)
        main_toolbar.addAction(self.open_action)
        main_toolbar.addAction(self.open_example_action)
        main_toolbar.addAction(self.save_action)
        main_toolbar.addSeparator()
        main_toolbar.addAction(self.solve_action)
        main_toolbar.addAction(self.export_plot_action)

    def _create_status_bar(self):
        """Create status bar."""
        self.status_label = QLabel("Ready")
        self.statusBar().addPermanentWidget(self.status_label)
        self.statusBar().showMessage("Ready", 3000)

    def _create_docks(self):
        """Create docked panels."""
        # Model tree dock
        self.tree_dock = QDockWidget("Model Tree", self)
        self.tree_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.model_tree = QTreeWidget()
        self.model_tree.setHeaderLabel("Model")
        self._populate_tree()
        self.tree_dock.setWidget(self.model_tree)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.tree_dock)

        # Properties dock
        self.props_dock = QDockWidget("Properties", self)
        self.props_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.props_panel = QTextEdit()
        self.props_panel.setReadOnly(True)
        self.props_panel.setPlaceholderText("Select an item to view properties...")
        self.props_dock.setWidget(self.props_panel)
        self.addDockWidget(Qt.RightDockWidgetArea, self.props_dock)

        # Visualization controls dock
        self.viz_ctrl_dock = QDockWidget("Visualization Controls", self)
        self.viz_ctrl_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        ctrl_widget = QWidget()
        ctrl_layout = QFormLayout(ctrl_widget)

        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setDecimals(2)
        self.scale_spin.setRange(0.01, 1000000.0)
        self.scale_spin.setValue(1000.0)
        self.scale_spin.setSingleStep(100.0)
        self.scale_spin.valueChanged.connect(self._refresh_visualization)

        self.reaction_scale_spin = QDoubleSpinBox()
        self.reaction_scale_spin.setDecimals(4)
        self.reaction_scale_spin.setRange(0.0001, 10.0)
        self.reaction_scale_spin.setValue(0.05)
        self.reaction_scale_spin.setSingleStep(0.01)
        self.reaction_scale_spin.valueChanged.connect(self._refresh_visualization)

        self.chk_deformed = QCheckBox("Show deformed shape")
        self.chk_deformed.setChecked(True)
        self.chk_deformed.toggled.connect(self.show_deformed_action.setChecked)
        self.chk_deformed.toggled.connect(self._refresh_visualization)

        self.chk_forces = QCheckBox("Color by axial force")
        self.chk_forces.setChecked(True)
        self.chk_forces.toggled.connect(self.show_forces_action.setChecked)
        self.chk_forces.toggled.connect(self._refresh_visualization)

        self.chk_reactions = QCheckBox("Show reactions")
        self.chk_reactions.setChecked(True)
        self.chk_reactions.toggled.connect(self.show_reactions_action.setChecked)
        self.chk_reactions.toggled.connect(self._refresh_visualization)

        self.btn_redraw = QPushButton("Redraw")
        self.btn_redraw.clicked.connect(self._refresh_visualization)

        ctrl_layout.addRow("Deformation scale", self.scale_spin)
        ctrl_layout.addRow("Reaction scale", self.reaction_scale_spin)
        ctrl_layout.addRow(self.chk_deformed)
        ctrl_layout.addRow(self.chk_forces)
        ctrl_layout.addRow(self.chk_reactions)
        ctrl_layout.addRow(self.btn_redraw)

        self.viz_ctrl_dock.setWidget(ctrl_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, self.viz_ctrl_dock)

    def _create_central_widget(self):
        """Create central widget with interactive visualization and run log."""
        central = QWidget()
        layout = QVBoxLayout()

        self.figure = Figure(figsize=(8, 6), dpi=110)
        self.canvas = FigureCanvas(self.figure)
        self.nav_toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        self._draw_empty_visualization()

        self.log_panel = QTextEdit()
        self.log_panel.setReadOnly(True)
        self.log_panel.setPlaceholderText("Run logs will appear here...")
        self.log_panel.setMaximumHeight(220)

        layout.addWidget(self.nav_toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.log_panel)
        central.setLayout(layout)
        self.setCentralWidget(central)

    def _draw_empty_visualization(self):
        """Draw placeholder on the matplotlib canvas."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title("Visualization Preview")
        self.ax.text(
            0.5,
            0.5,
            "1. Open a model file (txt)\n2. Click Solve (F5)\n3. Explore the plot with pan/zoom tools",
            ha="center",
            va="center",
            transform=self.ax.transAxes,
            color="#777777"
        )
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas.draw_idle()

    def _populate_tree(self):
        """Populate model tree with default items."""
        self.model_tree.clear()

        # Top-level items
        nodes_item = QTreeWidgetItem(self.model_tree, ["Nodes (0)"])
        elements_item = QTreeWidgetItem(self.model_tree, ["Elements (0)"])
        materials_item = QTreeWidgetItem(self.model_tree, ["Materials (0)"])
        sections_item = QTreeWidgetItem(self.model_tree, ["Sections (0)"])
        loads_item = QTreeWidgetItem(self.model_tree, ["Loads (0)"])
        constraints_item = QTreeWidgetItem(self.model_tree, ["Constraints (0)"])

        self.model_tree.expandAll()

    def _update_window_title(self):
        """Update window title based on current file."""
        title = "FEM_2d - 2D Finite Element Analysis"
        if self.current_file:
            title += f" - {os.path.basename(self.current_file)}"
        if self.is_modified:
            title += " *"
        self.setWindowTitle(title)

    def _update_actions(self):
        """Update action states based on current state."""
        has_model = bool(self.current_file)
        self.save_action.setEnabled(has_model and self.is_modified)
        self.save_as_action.setEnabled(has_model)
        self.solve_action.setEnabled(has_model)
        self.export_plot_action.setEnabled(self.results_data is not None)

    def _append_log(self, text):
        """Append one log line to the GUI log panel."""
        self.log_panel.append(text)

    def _refresh_visualization(self, *_args):
        """Re-render plot based on current controls and cached result data."""
        self.chk_deformed.blockSignals(True)
        self.chk_forces.blockSignals(True)
        self.chk_reactions.blockSignals(True)
        self.chk_deformed.setChecked(self.show_deformed_action.isChecked())
        self.chk_forces.setChecked(self.show_forces_action.isChecked())
        self.chk_reactions.setChecked(self.show_reactions_action.isChecked())
        self.chk_deformed.blockSignals(False)
        self.chk_forces.blockSignals(False)
        self.chk_reactions.blockSignals(False)

        if not self.results_data:
            self._draw_empty_visualization()
            return

        self._render_visualization(self.results_data)

    def _parse_results(self, results_path):
        """Parse results file produced by fem_run."""
        if not os.path.exists(results_path):
            raise FileNotFoundError(f"Results file not found: {results_path}")

        with open(results_path, "r", encoding="utf-8") as f:
            lines = [ln.rstrip("\n") for ln in f.readlines()]

        n_nodes = 0
        nodes = []
        elements = []
        displ = {}
        end_forces = {}
        reactions = {}

        idx = 0
        while idx < len(lines):
            line = lines[idx].strip()
            if line.startswith("Total nodes:"):
                m = re.search(r"Total nodes:\\s+(\\d+)", line)
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
                    nums = re.findall(r"[-+]?\\d*\\.?\\d+(?:e[-+]?\\d+)?", t, flags=re.IGNORECASE)
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
                    nums = re.findall(r"[-+]?\\d*\\.?\\d+(?:e[-+]?\\d+)?", t, flags=re.IGNORECASE)
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
                    nums = re.findall(r"[-+]?\\d*\\.?\\d+(?:e[-+]?\\d+)?", t, flags=re.IGNORECASE)
                    if len(nums) >= 4:
                        node = int(float(nums[0]))
                        rx = float(nums[1])
                        ry = float(nums[2])
                        rz = float(nums[3])
                        reactions[node] = (rx, ry, rz)
                    idx += 1
                continue

            idx += 1

        return {
            "nodes": nodes,
            "elements": elements,
            "displ": displ,
            "end_forces": end_forces,
            "reactions": reactions,
        }

    def _render_visualization(self, result):
        """Render structure and deformed plot on the canvas."""
        nodes = result["nodes"]
        elements = result["elements"]
        displ = result["displ"]
        end_forces = result["end_forces"]
        reactions = result["reactions"]

        if not nodes or not elements:
            self._draw_empty_visualization()
            self._append_log("[WARN] No node/element data to visualize")
            return

        scale = self.scale_spin.value()
        show_deformed = self.show_deformed_action.isChecked()
        show_forces = self.show_forces_action.isChecked()
        show_reactions = self.show_reactions_action.isChecked()
        reaction_scale = self.reaction_scale_spin.value()

        self.figure.clear()
        self.ax = self.figure.add_subplot(111)

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
            ux, uy, _ = displ.get(i, (0.0, 0.0, 0.0))
            dx = scale * ux
            dy = scale * uy
            if (not math.isfinite(dx)) or abs(dx) > cap:
                dx = 0.0
            if (not math.isfinite(dy)) or abs(dy) > cap:
                dy = 0.0
            xs1.append(nodes[i][1] + dx)
            ys1.append(nodes[i][2] + dy)

        axial_vals = []
        for idx_e, (_, i0, i1, _, _) in enumerate(elements):
            ef = end_forces.get(idx_e)
            if ef is None:
                axial_vals.append(0.0)
                continue
            x0, y0 = xs0[i0], ys0[i0]
            x1, y1 = xs0[i1], ys0[i1]
            dx0 = x1 - x0
            dy0 = y1 - y0
            length = math.hypot(dx0, dy0)
            if length <= 0:
                axial_vals.append(0.0)
                continue
            ux = dx0 / length
            uy = dy0 / length
            fx_i, fy_i = ef[0], ef[1]
            fx_j, fy_j = ef[3], ef[4]
            n_i = fx_i * ux + fy_i * uy
            n_j = fx_j * ux + fy_j * uy
            axial_vals.append(abs(0.5 * (n_i + n_j)))

        vmin = min(axial_vals) if axial_vals else 0.0
        vmax = max(axial_vals) if axial_vals else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colormaps["plasma"]

        for _, i0, i1, _, _ in elements:
            self.ax.plot([xs0[i0], xs0[i1]], [ys0[i0], ys0[i1]], color="#8b8f98", linewidth=1.8, alpha=0.9)
        self.ax.scatter(xs0, ys0, s=24, color="#59606d", label="original")

        if show_deformed:
            for idx_e, (_, i0, i1, _, _) in enumerate(elements):
                color = "#d62728"
                if show_forces and axial_vals:
                    color = cmap(norm(axial_vals[idx_e]))
                self.ax.plot([xs1[i0], xs1[i1]], [ys1[i0], ys1[i1]], color=color, linewidth=2.1)
            self.ax.scatter(xs1, ys1, s=26, color="#d62728", label="deformed")

        if show_reactions and reactions:
            arrow_scale = reaction_scale * base
            for node, (rx, ry, _) in reactions.items():
                x = xs0[node]
                y = ys0[node]
                self.ax.quiver(
                    [x],
                    [y],
                    [rx * arrow_scale],
                    [ry * arrow_scale],
                    angles="xy",
                    scale_units="xy",
                    scale=1.0,
                    color="#1f77b4",
                    width=0.0035,
                )

        if show_forces and axial_vals:
            sm = colormaps["plasma"]
            mappable = None
            try:
                from matplotlib.cm import ScalarMappable
                mappable = ScalarMappable(norm=norm, cmap=sm)
                mappable.set_array([])
            except Exception:
                mappable = None
            if mappable is not None:
                cbar = self.figure.colorbar(mappable, ax=self.ax, fraction=0.045, pad=0.03)
                cbar.set_label("|Axial Force|")

        self.ax.set_aspect("equal", adjustable="box")
        self.ax.grid(True, linestyle="--", alpha=0.25)
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_title(f"Structure / Deformed Visualization (scale={scale:g})")
        self.ax.legend(loc="upper right")
        self.figure.tight_layout()
        self.canvas.draw_idle()

    # Slots for actions
    @Slot()
    def new_model(self):
        """Create a new model."""
        # TODO: Check for unsaved changes
        self.model = None
        self.current_file = None
        self.is_modified = False
        self.results_data = None
        self.log_panel.clear()
        self._draw_empty_visualization()
        self._populate_tree()
        self._update_window_title()
        self._update_actions()
        self.statusBar().showMessage("New model created", 3000)

    @Slot()
    def open_default_example(self):
        """Open default example file from project root."""
        example = os.path.join(self.project_root, "test05.txt")
        if not os.path.exists(example):
            QMessageBox.warning(self, "Example Missing", f"Example file not found:\n{example}")
            return
        self.current_file = example
        self.model = object()
        self.is_modified = False
        self._update_window_title()
        self._update_actions()
        self.statusBar().showMessage("Opened default example: test05.txt", 3000)
        self._append_log(f"[INFO] Opened example: {example}")

    @Slot()
    def open_model(self):
        """Open an existing model."""
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open Model",
            "",
            "FEM Model Files (*.txt *.json);;All Files (*)"
        )

        if filename:
            self.current_file = filename
            self.model = object()
            self.is_modified = False
            self._update_window_title()
            self._update_actions()
            self.statusBar().showMessage(f"Loaded: {os.path.basename(filename)}", 3000)
            self._append_log(f"[INFO] Loaded model: {filename}")

    @Slot()
    def save_model(self):
        """Save the current model."""
        if not self.current_file:
            return self.save_model_as()

        # TODO: Save model
        self.is_modified = False
        self._update_window_title()
        self._update_actions()
        self.statusBar().showMessage("Model saved", 3000)
        return True

    @Slot()
    def save_model_as(self):
        """Save the model with a new name."""
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Model As",
            "",
            "JSON Files (*.json);;Text Files (*.txt);;All Files (*)"
        )

        if filename:
            self.current_file = filename
            return self.save_model()

        return False

    @Slot()
    def solve_model(self):
        """Solve model by calling release executable and render interactive plot."""
        if not self.current_file:
            QMessageBox.warning(self, "No Model", "Please open a model file first.")
            return

        exe_path = os.path.join(self.project_root, "build", "bin", "Release", "fem_run.exe")
        results_path = self.last_results_path

        if not os.path.exists(exe_path):
            QMessageBox.warning(
                self,
                "Executable Missing",
                "Release executable not found. Please build first:\n"
                "cmake --build build --config Release"
            )
            return

        os.makedirs(os.path.dirname(results_path), exist_ok=True)

        self.statusBar().showMessage("Solving...", 0)
        self._append_log("[INFO] Running solver...")

        try:
            solve_cmd = [exe_path, "--input", self.current_file, "--output", results_path]
            solve_proc = subprocess.run(
                solve_cmd,
                cwd=self.project_root,
                capture_output=True,
                text=True,
                check=False
            )
            if solve_proc.stdout:
                self._append_log(solve_proc.stdout)
            if solve_proc.stderr:
                self._append_log(solve_proc.stderr)
            if solve_proc.returncode != 0:
                self.statusBar().showMessage("Solve failed", 5000)
                QMessageBox.critical(self, "Solve Failed", "Solver execution failed. See logs for details.")
                return

            self.results_data = self._parse_results(results_path)
            self._render_visualization(self.results_data)
            self.figure.savefig(self.last_plot_path, dpi=150, bbox_inches="tight")
            self.statusBar().showMessage("Solution complete", 3000)
            self._append_log(f"[INFO] Results: {results_path}")
            self._append_log(f"[INFO] Plot: {self.last_plot_path}")
            self._update_actions()

        except Exception as exc:
            self.statusBar().showMessage("Solve failed", 5000)
            QMessageBox.critical(self, "Error", f"Unexpected error:\n{exc}")

    @Slot()
    def export_plot(self):
        """Export current visualization image to a user-selected path."""
        if self.results_data is None:
            QMessageBox.information(self, "No Plot", "Please run Solve first.")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Visualization",
            os.path.join(self.project_root, "build", "plot_export.png"),
            "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg);;SVG Files (*.svg);;All Files (*)"
        )
        if not filename:
            return

        try:
            self.figure.savefig(filename, dpi=200, bbox_inches="tight")
            self.statusBar().showMessage(f"Exported plot: {os.path.basename(filename)}", 4000)
            self._append_log(f"[INFO] Exported plot: {filename}")
        except Exception as exc:
            QMessageBox.critical(self, "Export Failed", f"Failed to export plot:\n{exc}")

    @Slot()
    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About FEM_2d",
            "<h3>FEM_2d</h3>"
            "<p>Version 1.0.0</p>"
            "<p>A 2D Finite Element Analysis tool for engineering education and practice.</p>"
            "<p>Supports truss and frame analysis with modern visualization.</p>"
            "<p><b>Features:</b></p>"
            "<ul>"
            "<li>Interactive model building</li>"
            "<li>Real-time visualization</li>"
            "<li>Internal force diagrams</li>"
            "<li>Educational examples</li>"
            "</ul>"
        )

    def closeEvent(self, event):
        """Handle window close event."""
        if self.is_modified:
            reply = QMessageBox.question(
                self,
                "Unsaved Changes",
                "The model has unsaved changes. Do you want to save before closing?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel
            )

            if reply == QMessageBox.Save:
                if not self.save_model():
                    event.ignore()
                    return
            elif reply == QMessageBox.Cancel:
                event.ignore()
                return

        event.accept()
