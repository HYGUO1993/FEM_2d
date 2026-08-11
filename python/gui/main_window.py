"""
Main Window for FEM_2d GUI Application

Modern layout: left side-panel (tabs: Model Tree | Editor | Controls) +
central canvas area + bottom log panel — all connected via QSplitter so
every region is freely resizable. Uses qtawesome for vector icons and a
QPropertyAnimation spinner during solve.
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QMenuBar, QMenu, QToolBar, QStatusBar, QLabel,
    QFileDialog, QMessageBox, QSplitter, QPushButton,
    QCheckBox, QDoubleSpinBox, QFormLayout, QPlainTextEdit,
    QTabWidget, QTreeWidget, QTreeWidgetItem, QTextEdit,
    QSizePolicy, QFrame, QGraphicsDropShadowEffect,
)
from PyQt6.QtCore import (
    Qt, pyqtSignal as Signal, pyqtSlot as Slot,
    QPropertyAnimation, QEasingCurve, QSize, QThread,
    pyqtProperty,
)
from PyQt6.QtGui import QAction, QKeySequence, QFont, QColor, QTransform

import os
import subprocess
import tempfile

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

try:
    import qtawesome as qta
    _HAS_QTA = True
except ImportError:
    _HAS_QTA = False

from .theme import setup_matplotlib_dark_theme, PALETTE
from .fem_parser import parse_input_file, parse_results_file, FEMModelData
from .visualization import FEMVisualizer


# ---------------------------------------------------------------------------
# Helper: build a qtawesome icon or fall back to empty icon
# ---------------------------------------------------------------------------
def _icon(name: str, color: str = PALETTE["text"]):
    if _HAS_QTA:
        return qta.icon(name, color=color)
    from PyQt6.QtGui import QIcon
    return QIcon()


def _accent_button(text: str, icon_name: str = None) -> QPushButton:
    """QPushButton styled as a primary accent action."""
    btn = QPushButton(text)
    btn.setProperty("accent", "true")
    if icon_name:
        btn.setIcon(_icon(icon_name, color=PALETTE["crust"]))
    btn.setMinimumHeight(32)
    return btn


def _secondary_button(text: str, icon_name: str = None) -> QPushButton:
    btn = QPushButton(text)
    if icon_name:
        btn.setIcon(_icon(icon_name))
    btn.setMinimumHeight(28)
    return btn


# ---------------------------------------------------------------------------
# Solve worker thread
# ---------------------------------------------------------------------------
class SolveWorker(QThread):
    finished = Signal(bool, str)
    log_line = Signal(str)

    def __init__(self, exe_path, input_path, output_path, work_dir):
        super().__init__()
        self._exe = exe_path
        self._input = input_path
        self._output = output_path
        self._work_dir = work_dir

    def run(self):
        try:
            info = None
            if os.name == "nt":
                info = subprocess.STARTUPINFO()
                info.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            cmd = [self._exe, "--input", self._input, "--output", self._output, "--quiet"]
            proc = subprocess.run(
                cmd, cwd=self._work_dir,
                capture_output=True, text=True, startupinfo=info,
            )
            if proc.stdout:
                self.log_line.emit(proc.stdout.strip())
            if proc.stderr:
                self.log_line.emit("[ERR] " + proc.stderr.strip())
            if proc.returncode != 0:
                self.finished.emit(False, "求解器返回非零退出码，请查看日志。")
            else:
                self.finished.emit(True, "求解完成！")
        except Exception as exc:
            self.finished.emit(False, str(exc))


# ---------------------------------------------------------------------------
# Animated spinner label
# ---------------------------------------------------------------------------
class SpinnerLabel(QLabel):
    """Rotating icon used as a progress indicator during solve."""

    def __init__(self, parent=None):
        super().__init__(parent)
        if _HAS_QTA:
            self._base_pixmap = qta.icon(
                "fa5s.circle-notch", color=PALETTE["blue"]
            ).pixmap(QSize(18, 18))
        else:
            self._base_pixmap = None
        self._angle = 0.0
        self._anim = QPropertyAnimation(self, b"angle")
        self._anim.setStartValue(0.0)
        self._anim.setEndValue(360.0)
        self._anim.setDuration(900)
        self._anim.setLoopCount(-1)
        self._anim.setEasingCurve(QEasingCurve.Type.Linear)
        self.setFixedSize(20, 20)
        self.hide()

    def _get_angle(self):
        return self._angle

    def _set_angle(self, value):
        self._angle = value
        if self._base_pixmap:
            t = QTransform().rotate(value)
            self.setPixmap(
                self._base_pixmap.transformed(t, Qt.TransformationMode.SmoothTransformation)
            )

    angle = pyqtProperty(float, _get_angle, _set_angle)

    def start(self):
        self.show()
        self._anim.start()

    def stop(self):
        self._anim.stop()
        self.hide()


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------
class MainWindow(QMainWindow):
    """Main application window — VS Code-style splitter layout."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("FEM_2d — 2D Finite Element Analysis")
        self.setGeometry(100, 80, 1440, 900)
        self.setMinimumSize(960, 640)

        self.current_file = None
        self.model = None
        self.is_modified = False
        self._solve_worker = None

        import sys
        if getattr(sys, "frozen", False):
            self.app_dir = os.path.dirname(sys.executable)
        else:
            self.app_dir = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "..", "..")
            )
        self.work_dir = self.app_dir
        self.last_results_path = os.path.join(tempfile.gettempdir(), "FEM_2d_Results.dat")

        setup_matplotlib_dark_theme()

        self._create_actions()
        self._create_menus()
        self._create_toolbar()
        self._create_status_bar()
        self._create_layout()

        self._update_window_title()
        self._update_actions()

    # ── Actions ──────────────────────────────────────────────────────────
    def _create_actions(self):
        def act(label, shortcut=None, icon=None, checkable=False):
            a = QAction(label, self)
            if shortcut:
                a.setShortcut(shortcut)
            if icon:
                a.setIcon(_icon(icon))
            if checkable:
                a.setCheckable(True)
            return a

        self.new_action    = act("新建 (&N)", QKeySequence.StandardKey.New,  "fa5s.file")
        self.open_action   = act("打开... (&O)", QKeySequence.StandardKey.Open, "fa5s.folder-open")
        self.save_action   = act("保存 (&S)", QKeySequence.StandardKey.Save,  "fa5s.save")
        self.exit_action   = act("退出 (&X)", QKeySequence.StandardKey.Quit)
        self.solve_action  = act("求解运行 (F5)", Qt.Key.Key_F5, "fa5s.play")
        self.open_example_action = act("加载演示算例")
        self.export_plot_action  = act("导出图片...", None, "fa5s.image")
        self.about_action        = act("关于")

        self.show_deformed_action  = act("显示变形",  checkable=True)
        self.show_forces_action    = act("内力色图",  checkable=True)
        self.show_reactions_action = act("支座反力",  checkable=True)
        self.show_deformed_action.setChecked(True)
        self.show_forces_action.setChecked(True)
        self.show_reactions_action.setChecked(True)

        self.new_action.triggered.connect(self.new_model)
        self.open_action.triggered.connect(self.open_model)
        self.save_action.triggered.connect(self.save_model)
        self.exit_action.triggered.connect(self.close)
        self.solve_action.triggered.connect(self.solve_model)
        self.open_example_action.triggered.connect(self.open_default_example)
        self.export_plot_action.triggered.connect(self.export_plot)
        self.about_action.triggered.connect(self.show_about)
        self.show_deformed_action.triggered.connect(self._refresh_visualization)
        self.show_forces_action.triggered.connect(self._refresh_visualization)
        self.show_reactions_action.triggered.connect(self._refresh_visualization)

    # ── Menus ─────────────────────────────────────────────────────────────
    def _create_menus(self):
        mb = self.menuBar()

        fm = mb.addMenu("文件 (&F)")
        fm.addAction(self.new_action)
        fm.addAction(self.open_action)
        fm.addAction(self.open_example_action)
        fm.addSeparator()
        fm.addAction(self.save_action)
        fm.addSeparator()
        fm.addAction(self.exit_action)

        am = mb.addMenu("分析 (&A)")
        am.addAction(self.solve_action)
        am.addAction(self.export_plot_action)

        vm = mb.addMenu("视图 (&V)")
        vm.addAction(self.show_deformed_action)
        vm.addAction(self.show_forces_action)
        vm.addAction(self.show_reactions_action)

        hm = mb.addMenu("帮助 (&H)")
        hm.addAction(self.about_action)

    # ── Toolbar ───────────────────────────────────────────────────────────
    def _create_toolbar(self):
        tb = self.addToolBar("主工具栏")
        tb.setMovable(False)
        tb.setIconSize(QSize(18, 18))
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)

        tb.addAction(self.new_action)
        tb.addAction(self.open_action)
        tb.addAction(self.save_action)
        tb.addSeparator()
        tb.addAction(self.solve_action)
        tb.addAction(self.export_plot_action)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)
        self._spinner = SpinnerLabel()
        tb.addWidget(self._spinner)

    # ── Status bar ────────────────────────────────────────────────────────
    def _create_status_bar(self):
        self._status_file_label = QLabel("就绪")
        self._status_file_label.setProperty("subheading", "true")
        self.statusBar().addPermanentWidget(self._status_file_label)

    # ── Main layout ───────────────────────────────────────────────────────
    def _create_layout(self):
        root = QSplitter(Qt.Orientation.Horizontal)
        root.setHandleWidth(3)
        root.setChildrenCollapsible(False)

        self._side_tabs = QTabWidget()
        self._side_tabs.setTabPosition(QTabWidget.TabPosition.North)
        self._side_tabs.setMinimumWidth(220)
        self._side_tabs.setMaximumWidth(380)

        self._side_tabs.addTab(self._build_tree_panel(),     _icon("fa5s.sitemap",   PALETTE["blue"]),  "模型树")
        self._side_tabs.addTab(self._build_editor_panel(),   _icon("fa5s.edit",      PALETTE["mauve"]), "编辑")
        self._side_tabs.addTab(self._build_controls_panel(), _icon("fa5s.sliders-h", PALETTE["teal"]),  "控制")

        right_split = QSplitter(Qt.Orientation.Vertical)
        right_split.setHandleWidth(3)
        right_split.setChildrenCollapsible(False)
        right_split.addWidget(self._build_canvas_widget())
        right_split.addWidget(self._build_log_widget())
        right_split.setStretchFactor(0, 5)
        right_split.setStretchFactor(1, 1)
        right_split.setSizes([700, 130])

        root.addWidget(self._side_tabs)
        root.addWidget(right_split)
        root.setStretchFactor(0, 0)
        root.setStretchFactor(1, 1)
        root.setSizes([260, 1000])

        self.setCentralWidget(root)

    # -- Side panel tabs --------------------------------------------------
    def _build_tree_panel(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)

        hdr = QLabel("模型结构")
        hdr.setProperty("heading", "true")
        lay.addWidget(hdr)

        self.model_tree = QTreeWidget()
        self.model_tree.setHeaderLabel("组件")
        self.model_tree.setAlternatingRowColors(True)
        self._populate_tree()

        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(12)
        shadow.setOffset(0, 2)
        shadow.setColor(QColor(0, 0, 0, 80))
        self.model_tree.setGraphicsEffect(shadow)

        lay.addWidget(self.model_tree)
        return w

    def _build_editor_panel(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)

        hdr = QLabel("输入数据编辑器")
        hdr.setProperty("heading", "true")
        lay.addWidget(hdr)

        sub = QLabel("编辑后点击「解析应用」更新模型")
        sub.setProperty("subheading", "true")
        sub.setWordWrap(True)
        lay.addWidget(sub)

        self.text_editor = QPlainTextEdit()
        self.text_editor.setFont(QFont("Consolas", 10))
        self.text_editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.text_editor.textChanged.connect(self._mark_modified)
        lay.addWidget(self.text_editor, stretch=1)

        btn_apply = _accent_button("⚙  解析并应用", "fa5s.check")
        btn_apply.clicked.connect(self._apply_editor_text)
        lay.addWidget(btn_apply)
        return w

    def _build_controls_panel(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(8)

        hdr = QLabel("视图控制")
        hdr.setProperty("heading", "true")
        lay.addWidget(hdr)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        form.setSpacing(8)

        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.1, 1_000_000.0)
        self.scale_spin.setValue(100.0)
        self.scale_spin.setSingleStep(100.0)
        self.scale_spin.setDecimals(1)
        self.scale_spin.setSuffix("×")
        self.scale_spin.valueChanged.connect(self._refresh_visualization)
        form.addRow("变形放大", self.scale_spin)

        self.reaction_scale_spin = QDoubleSpinBox()
        self.reaction_scale_spin.setRange(0.001, 10.0)
        self.reaction_scale_spin.setValue(0.05)
        self.reaction_scale_spin.setSingleStep(0.01)
        self.reaction_scale_spin.setDecimals(3)
        self.reaction_scale_spin.valueChanged.connect(self._refresh_visualization)
        form.addRow("反力比例", self.reaction_scale_spin)

        lay.addLayout(form)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet(f"background-color: {PALETTE['surface0']}; max-height: 1px;")
        lay.addWidget(sep)

        self.chk_nodes     = QCheckBox("显示节点编号")
        self.chk_elems     = QCheckBox("显示单元编号")
        self.chk_deformed  = QCheckBox("显示位移变形")
        self.chk_forces    = QCheckBox("轴力色图")
        self.chk_bmd       = QCheckBox("弯矩图 (BMD)")
        self.chk_reactions = QCheckBox("显示支座反力")

        self.chk_deformed.setChecked(True)
        self.chk_forces.setChecked(True)
        self.chk_reactions.setChecked(True)

        for chk in (self.chk_nodes, self.chk_elems, self.chk_deformed,
                    self.chk_forces, self.chk_bmd, self.chk_reactions):
            lay.addWidget(chk)
            chk.toggled.connect(self._refresh_visualization)

        self.chk_deformed.toggled.connect(self.show_deformed_action.setChecked)
        self.show_deformed_action.toggled.connect(self.chk_deformed.setChecked)
        self.chk_forces.toggled.connect(self.show_forces_action.setChecked)
        self.show_forces_action.toggled.connect(self.chk_forces.setChecked)
        self.chk_reactions.toggled.connect(self.show_reactions_action.setChecked)
        self.show_reactions_action.toggled.connect(self.chk_reactions.setChecked)

        sep2 = QFrame()
        sep2.setFrameShape(QFrame.Shape.HLine)
        sep2.setStyleSheet(f"background-color: {PALETTE['surface0']}; max-height: 1px;")
        lay.addWidget(sep2)

        btn_redraw = _secondary_button("↺  刷新画板", "fa5s.sync-alt")
        btn_redraw.clicked.connect(self._refresh_visualization)
        lay.addWidget(btn_redraw)

        btn_export = _secondary_button("📸  导出图片", "fa5s.image")
        btn_export.clicked.connect(self.export_plot)
        lay.addWidget(btn_export)

        btn_report = _secondary_button("📄  打开计算书", "fa5s.file-alt")
        btn_report.clicked.connect(self.open_text_report)
        lay.addWidget(btn_report)

        lay.addStretch()
        return w

    def _build_canvas_widget(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 0)
        lay.setSpacing(0)

        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        self.figure = Figure(figsize=(8, 6), dpi=110)
        self.figure.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.10)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._nav_toolbar = NavigationToolbar2QT(self.canvas, None)
        self._nav_toolbar.setIconSize(QSize(16, 16))

        self.visualizer = FEMVisualizer(self.figure, self.canvas)
        self.visualizer.draw_empty()

        lay.addWidget(self._nav_toolbar)
        lay.addWidget(self.canvas, stretch=1)
        return w

    def _build_log_widget(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 0, 4, 4)
        lay.setSpacing(2)

        hdr_row = QHBoxLayout()
        lbl = QLabel("  日志输出")
        lbl.setProperty("subheading", "true")
        hdr_row.addWidget(lbl)
        hdr_row.addStretch()
        btn_clear = QPushButton("清空")
        btn_clear.setFixedWidth(52)
        btn_clear.setFixedHeight(22)
        btn_clear.clicked.connect(lambda: self.log_panel.clear())
        hdr_row.addWidget(btn_clear)
        lay.addLayout(hdr_row)

        self.log_panel = QTextEdit()
        self.log_panel.setReadOnly(True)
        self.log_panel.setPlaceholderText("系统日志输出…")
        self.log_panel.setFont(QFont("Consolas", 9))
        self.log_panel.setMaximumHeight(160)
        lay.addWidget(self.log_panel)
        return w

    # ── Model tree ────────────────────────────────────────────────────────
    def _populate_tree(self):
        self.model_tree.clear()
        if not self.model:
            for label in ["Nodes (0)", "Elements (0)", "Materials (0)",
                          "Sections (0)", "Loads (0)", "Constraints (0)"]:
                QTreeWidgetItem(self.model_tree, [label])
            return

        def add_group(label, children):
            item = QTreeWidgetItem(self.model_tree, [label])
            for c in children:
                QTreeWidgetItem(item, [c])
            return item

        nodes_item = add_group(
            f"Nodes ({len(self.model.nodes)})",
            [f"Node {i}: {'Frame' if t==2 else 'Truss'}  ({x}, {y})"
             for i, (t, x, y) in enumerate(self.model.nodes)],
        )
        elems_item = add_group(
            f"Elements ({len(self.model.elements)})",
            [f"Elem {i}: N{i0}→N{i1}  (Sec:{sec}, Mat:{mat})"
             for i, (_, i0, i1, sec, mat) in enumerate(self.model.elements)],
        )
        add_group(
            f"Materials ({len(self.model.materials)})",
            [f"Mat {i}: E={e:g}, v={mu:g}" for i, (e, mu, _) in enumerate(self.model.materials)],
        )
        add_group(
            f"Sections ({len(self.model.sections)})",
            [f"Sec {i}: A={A:g}, Iz={Iz:g}" for i, (A, Iz, _) in enumerate(self.model.sections)],
        )
        add_group(
            f"Constraints ({len(self.model.constraints)})",
            [f"Node {n}: Fixed({''.join(k for k,v in zip(['X','Y','Rz'],[cx,cy,cr]) if v<0)})"
             for n, (cx, cy, cr) in self.model.constraints.items()],
        )
        add_group(
            f"Loads ({len(self.model.loads)})",
            [f"{'Node' if lt==1 else 'Dist.'} load — {['X','Y','Rz'][ld] if ld<3 else ld}={val:g}"
             for lt, ld, val, *_ in self.model.loads],
        )

        for item in (nodes_item, elems_item):
            self.model_tree.expandItem(item)

    # ── State helpers ─────────────────────────────────────────────────────
    def _update_window_title(self):
        title = "FEM_2d"
        if self.current_file:
            title += f"  —  {os.path.basename(self.current_file)}"
        if self.is_modified:
            title += "  ●"
        self.setWindowTitle(title)

    def _update_actions(self):
        has = self.model is not None
        self.save_action.setEnabled(has and self.is_modified)
        self.solve_action.setEnabled(has)
        self.export_plot_action.setEnabled(has and getattr(self.model, "is_solved", False))

    def _append_log(self, text):
        self.log_panel.append(text)

    def _mark_modified(self):
        if not self.is_modified:
            self.is_modified = True
            self._update_window_title()
            self._update_actions()

    @Slot()
    def _refresh_visualization(self):
        if not self.model:
            self.visualizer.draw_empty()
            return
        self.visualizer.render(
            self.model,
            scale=self.scale_spin.value(),
            reaction_scale=self.reaction_scale_spin.value(),
            show_deformed=self.chk_deformed.isChecked(),
            show_forces=self.chk_forces.isChecked(),
            show_reactions=self.chk_reactions.isChecked(),
            show_loads=True,
            show_nodes=self.chk_nodes.isChecked(),
            show_elems=self.chk_elems.isChecked(),
            show_bmd=self.chk_bmd.isChecked(),
        )

    # ── File operations ───────────────────────────────────────────────────
    def _load_model_file(self, filepath):
        try:
            self.model = parse_input_file(filepath)
            self.current_file = filepath
            self.is_modified = False

            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
            self.text_editor.blockSignals(True)
            self.text_editor.setPlainText(content)
            self.text_editor.blockSignals(False)

            self._populate_tree()
            self._update_window_title()
            self._update_actions()
            self._refresh_visualization()
            self._side_tabs.setCurrentIndex(0)

            msg = f"已加载: {os.path.basename(filepath)}"
            self.statusBar().showMessage(msg, 4000)
            self._status_file_label.setText(os.path.basename(filepath))
            self._append_log(f"[INFO] 成功加载模型: {filepath}")
        except Exception as exc:
            QMessageBox.critical(self, "加载失败", f"无法解析文件:\n{exc}")
            self._append_log(f"[ERROR] 读取失败: {exc}")

    @Slot()
    def _apply_editor_text(self):
        if not self.current_file:
            self.current_file = os.path.join(tempfile.gettempdir(), "fem2d_temp_input.txt")
        try:
            with open(self.current_file, "w", encoding="utf-8") as f:
                f.write(self.text_editor.toPlainText())
            self.model = parse_input_file(self.current_file)
            self.is_modified = False
            self._populate_tree()
            self._refresh_visualization()
            self._update_window_title()
            self._update_actions()
            self._append_log("[INFO] 已从编辑器文本解析模型")
            self.statusBar().showMessage("模型解析成功", 3000)
        except Exception as exc:
            QMessageBox.critical(self, "解析失败", f"文本格式有误:\n{exc}")

    @Slot()
    def new_model(self):
        self.model = None
        self.current_file = None
        self.is_modified = False
        self.text_editor.blockSignals(True)
        self.text_editor.clear()
        self.text_editor.blockSignals(False)
        self.log_panel.clear()
        self.visualizer.draw_empty()
        self._populate_tree()
        self._update_window_title()
        self._update_actions()
        self._status_file_label.setText("就绪")
        self.statusBar().showMessage("新建模型", 3000)

    @Slot()
    def open_default_example(self):
        example = os.path.join(self.app_dir, "example_continuous_beam.txt")
        if not os.path.exists(example):
            QMessageBox.warning(self, "文件丢失", f"找不到演示算例:\n{example}")
            return
        self._load_model_file(example)

    @Slot()
    def open_model(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "打开文件", self.work_dir, "文本文件 (*.txt);;所有文件 (*)"
        )
        if filename:
            self._load_model_file(filename)

    @Slot()
    def save_model(self):
        if not self.current_file:
            filename, _ = QFileDialog.getSaveFileName(
                self, "保存文件", self.work_dir, "文本文件 (*.txt);;所有文件 (*)"
            )
            if not filename:
                return False
            self.current_file = filename
        try:
            with open(self.current_file, "w", encoding="utf-8") as f:
                f.write(self.text_editor.toPlainText())
            self.is_modified = False
            self._update_window_title()
            self._update_actions()
            self.statusBar().showMessage("已保存", 3000)
            return True
        except Exception as exc:
            QMessageBox.critical(self, "保存失败", f"无法写入文件:\n{exc}")
            return False

    # ── Solver ────────────────────────────────────────────────────────────
    def _find_solver(self):
        import sys
        if getattr(sys, "frozen", False):
            return os.path.join(self.app_dir, "bin", "fem_run.exe")
        for candidate in [
            os.path.join(self.app_dir, "build", "bin", "Release", "fem_run.exe"),
            os.path.join(self.app_dir, "build", "bin", "Debug",   "fem_run.exe"),
            os.path.join(self.app_dir, "build", "bin", "Release", "fem_run"),
            os.path.join(self.app_dir, "build", "bin", "Debug",   "fem_run"),
        ]:
            if os.path.exists(candidate):
                return candidate
        return None

    @Slot()
    def solve_model(self):
        if not self.model:
            QMessageBox.warning(self, "无模型", "请先打开或创建模型。")
            return
        if self.is_modified:
            self.save_model()

        exe = self._find_solver()
        if not exe:
            QMessageBox.critical(
                self, "找不到求解器",
                "未找到 fem_run 可执行文件，请先编译项目。\n"
                "（Release 或 Debug 版本均可）",
            )
            return

        os.makedirs(os.path.dirname(self.last_results_path), exist_ok=True)
        self.statusBar().showMessage("求解中…", 0)
        self._append_log(f"[INFO] 启动求解器: {exe}")
        self._spinner.start()
        self.solve_action.setEnabled(False)

        self._solve_worker = SolveWorker(exe, self.current_file, self.last_results_path, self.work_dir)
        self._solve_worker.log_line.connect(self._append_log)
        self._solve_worker.finished.connect(self._on_solve_finished)
        self._solve_worker.start()

    @Slot(bool, str)
    def _on_solve_finished(self, success, message):
        self._spinner.stop()
        self.solve_action.setEnabled(True)

        if not success:
            QMessageBox.critical(self, "求解失败", message)
            self.statusBar().showMessage("求解失败", 5000)
            return

        try:
            self._append_log("[INFO] 解析结果文件…")
            self.model = parse_results_file(self.last_results_path, self.model)
            self._refresh_visualization()
            self._update_actions()
            self.statusBar().showMessage(message, 4000)
            self._append_log(f"[INFO] {message}")
        except Exception as exc:
            QMessageBox.critical(self, "结果解析失败", str(exc))
            self._append_log(f"[ERROR] 结果解析失败: {exc}")

    # ── Export / report ───────────────────────────────────────────────────
    @Slot()
    def export_plot(self):
        if not self.model or not self.model.is_solved:
            QMessageBox.information(self, "无结果", "请先求解模型再导出。")
            return
        default = os.path.join(os.path.expanduser("~"), "Desktop", "fem2d_result.png")
        filename, _ = QFileDialog.getSaveFileName(
            self, "导出图片", default, "PNG 图片 (*.png);;所有文件 (*)"
        )
        if filename:
            try:
                self.figure.savefig(filename, dpi=200, bbox_inches="tight")
                self.statusBar().showMessage(f"图片已导出: {os.path.basename(filename)}", 4000)
                self._append_log(f"[INFO] 已导出: {filename}")
            except Exception as exc:
                QMessageBox.critical(self, "导出失败", str(exc))

    @Slot()
    def open_text_report(self):
        if not os.path.exists(self.last_results_path):
            QMessageBox.information(self, "无计算书", "当前会话尚未完成求解，无结果文件。")
            return
        try:
            if os.name == "nt":
                os.startfile(self.last_results_path)
            else:
                import sys
                opener = "open" if sys.platform == "darwin" else "xdg-open"
                subprocess.run([opener, self.last_results_path])
        except Exception as exc:
            QMessageBox.critical(self, "打开失败", f"无法打开文件:\n{exc}")

    @Slot()
    def show_about(self):
        QMessageBox.about(
            self, "关于 FEM_2d",
            "<h3>FEM_2d</h3>"
            "<p>版本 2.1 — 现代化 UI 重构</p>"
            "<p>二维杆系与框架有限元分析程序</p>"
            "<p><small>计算内核: C++ · 界面: PyQt6 · 绘图: Matplotlib</small></p>",
        )

    # ── Close event ───────────────────────────────────────────────────────
    def closeEvent(self, event):
        if self.is_modified:
            reply = QMessageBox.question(
                self, "未保存修改",
                "当前有未保存的修改，退出前是否保存？",
                QMessageBox.StandardButton.Save |
                QMessageBox.StandardButton.Discard |
                QMessageBox.StandardButton.Cancel,
            )
            if reply == QMessageBox.StandardButton.Save:
                if not self.save_model():
                    event.ignore()
                    return
            elif reply == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
        if self._solve_worker and self._solve_worker.isRunning():
            self._solve_worker.terminate()
            self._solve_worker.wait(2000)
        event.accept()

