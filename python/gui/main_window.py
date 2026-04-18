"""
Main Window for FEM_2d GUI Application

This module implements the main application window with menu bar, toolbar,
status bar, and docked panels for model tree and properties.
"""

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QMenuBar, QMenu, QToolBar, QStatusBar, QDockWidget,
    QTreeWidget, QTreeWidgetItem, QTextEdit, QLabel,
    QFileDialog, QMessageBox, QSplitter, QPushButton,
    QCheckBox, QDoubleSpinBox, QFormLayout, QPlainTextEdit,
    QTabWidget
)
from PyQt6.QtCore import Qt, pyqtSignal as Signal, pyqtSlot as Slot
from PyQt6.QtGui import QAction, QIcon, QKeySequence, QFont

import os
import subprocess
import glob

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from .theme import setup_matplotlib_dark_theme
from .fem_parser import parse_input_file, parse_results_file, FEMModelData
from .visualization import FEMVisualizer

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
        import sys, tempfile
        if getattr(sys, 'frozen', False):
            self.app_dir = os.path.dirname(sys.executable)
            self.work_dir = os.getcwd()
        else:
            self.app_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            self.work_dir = self.app_dir

        self.last_results_path = os.path.join(tempfile.gettempdir(), "FEM_2d_Results.dat")
        self.last_plot_path = os.path.join(os.path.expanduser("~"), "Desktop", "plot_gui.png")

        # Setup matplotlib theme
        setup_matplotlib_dark_theme()

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
        self.new_action = QAction("新建 (&N)", self)
        self.new_action.setShortcut(QKeySequence.StandardKey.New)
        self.new_action.triggered.connect(self.new_model)

        self.open_action = QAction("打开... (&O)", self)
        self.open_action.setShortcut(QKeySequence.StandardKey.Open)
        self.open_action.triggered.connect(self.open_model)

        self.save_action = QAction("保存 (&S)", self)
        self.save_action.setShortcut(QKeySequence.StandardKey.Save)
        self.save_action.triggered.connect(self.save_model)

        self.exit_action = QAction("退出 (&X)", self)
        self.exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        self.exit_action.triggered.connect(self.close)

        # Analysis actions
        self.solve_action = QAction("求解运行 (F5)", self)
        self.solve_action.setShortcut(Qt.Key.Key_F5)
        self.solve_action.triggered.connect(self.solve_model)

        self.open_example_action = QAction("加载测试例子 (test05.txt)", self)
        self.open_example_action.triggered.connect(self.open_default_example)

        # View actions
        self.show_deformed_action = QAction("显示变形", self)
        self.show_deformed_action.setCheckable(True)
        self.show_deformed_action.setChecked(True)
        self.show_deformed_action.triggered.connect(self._refresh_visualization)

        self.show_forces_action = QAction("内力云图", self)
        self.show_forces_action.setCheckable(True)
        self.show_forces_action.setChecked(True)
        self.show_forces_action.triggered.connect(self._refresh_visualization)

        self.show_reactions_action = QAction("显示支座反力", self)
        self.show_reactions_action.setCheckable(True)
        self.show_reactions_action.setChecked(True)
        self.show_reactions_action.triggered.connect(self._refresh_visualization)

        self.export_plot_action = QAction("导出图片...", self)
        self.export_plot_action.triggered.connect(self.export_plot)

        self.about_action = QAction("关于", self)
        self.about_action.triggered.connect(self.show_about)

    def _create_menus(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("文件 (&F)")
        file_menu.addAction(self.new_action)
        file_menu.addAction(self.open_action)
        file_menu.addAction(self.open_example_action)
        file_menu.addSeparator()
        file_menu.addAction(self.save_action)
        file_menu.addSeparator()
        file_menu.addAction(self.exit_action)

        analysis_menu = menubar.addMenu("分析 (&A)")
        analysis_menu.addAction(self.solve_action)
        analysis_menu.addAction(self.export_plot_action)

        view_menu = menubar.addMenu("视图 (&V)")
        view_menu.addAction(self.show_deformed_action)
        view_menu.addAction(self.show_forces_action)
        view_menu.addAction(self.show_reactions_action)

        help_menu = menubar.addMenu("帮助 (&H)")
        help_menu.addAction(self.about_action)

    def _create_toolbars(self):
        main_toolbar = self.addToolBar("Main")
        main_toolbar.setMovable(False)
        main_toolbar.addAction(self.new_action)
        main_toolbar.addAction(self.open_action)
        main_toolbar.addAction(self.save_action)
        main_toolbar.addSeparator()
        main_toolbar.addAction(self.solve_action)
        main_toolbar.addAction(self.export_plot_action)

    def _create_status_bar(self):
        self.status_label = QLabel("就绪")
        self.statusBar().addPermanentWidget(self.status_label)
        self.statusBar().showMessage("就绪", 3000)

    def _create_docks(self):
        # Model tree dock
        self.tree_dock = QDockWidget("模型结构", self)
        self.tree_dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea)
        self.model_tree = QTreeWidget()
        self.model_tree.setHeaderLabel("组件")
        self._populate_tree()
        self.tree_dock.setWidget(self.model_tree)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.tree_dock)

        # File Editor Dock
        self.editor_dock = QDockWidget("输入数据编辑", self)
        self.editor_dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea | Qt.DockWidgetArea.RightDockWidgetArea)
        editor_widget = QWidget()
        el = QVBoxLayout(editor_widget)
        self.text_editor = QPlainTextEdit()
        font = QFont("Consolas", 10)
        self.text_editor.setFont(font)
        self.text_editor.textChanged.connect(self._mark_modified)
        
        btn_apply = QPushButton("解析并应用")
        btn_apply.clicked.connect(self._apply_editor_text)
        
        el.addWidget(self.text_editor)
        el.addWidget(btn_apply)
        self.editor_dock.setWidget(editor_widget)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.editor_dock)
        
        # Viz controls
        self.viz_ctrl_dock = QDockWidget("视图控制", self)
        self.viz_ctrl_dock.setAllowedAreas(Qt.DockWidgetArea.RightDockWidgetArea)
        ctrl_widget = QWidget()
        ctrl_layout = QFormLayout(ctrl_widget)

        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setDecimals(1)
        self.scale_spin.setRange(0.1, 1000000.0)
        self.scale_spin.setValue(100.0)
        self.scale_spin.setSingleStep(100.0)
        self.scale_spin.valueChanged.connect(self._refresh_visualization)

        self.reaction_scale_spin = QDoubleSpinBox()
        self.reaction_scale_spin.setDecimals(3)
        self.reaction_scale_spin.setRange(0.001, 10.0)
        self.reaction_scale_spin.setValue(0.05)
        self.reaction_scale_spin.setSingleStep(0.01)
        self.reaction_scale_spin.valueChanged.connect(self._refresh_visualization)

        self.chk_nodes = QCheckBox("显示节点编号")
        self.chk_nodes.setChecked(False)
        self.chk_nodes.toggled.connect(self._refresh_visualization)

        self.chk_elems = QCheckBox("显示单元编号")
        self.chk_elems.setChecked(False)
        self.chk_elems.toggled.connect(self._refresh_visualization)

        self.chk_deformed = QCheckBox("显示位移变形")
        self.chk_deformed.setChecked(True)
        self.chk_deformed.toggled.connect(self.show_deformed_action.setChecked)
        self.chk_deformed.toggled.connect(self._refresh_visualization)
        self.show_deformed_action.toggled.connect(self.chk_deformed.setChecked)

        self.chk_forces = QCheckBox("轴力颜色渲染")
        self.chk_forces.setChecked(True)
        self.chk_forces.toggled.connect(self.show_forces_action.setChecked)
        self.chk_forces.toggled.connect(self._refresh_visualization)
        self.show_forces_action.toggled.connect(self.chk_forces.setChecked)

        self.chk_bmd = QCheckBox("显示弯矩图 (BMD)")
        self.chk_bmd.setChecked(False)
        self.chk_bmd.toggled.connect(self._refresh_visualization)

        self.chk_reactions = QCheckBox("显示支座反力")
        self.chk_reactions.setChecked(True)
        self.chk_reactions.toggled.connect(self.show_reactions_action.setChecked)
        self.chk_reactions.toggled.connect(self._refresh_visualization)
        self.show_reactions_action.toggled.connect(self.chk_reactions.setChecked)

        self.btn_redraw = QPushButton("刷新画板")
        self.btn_redraw.clicked.connect(self._refresh_visualization)

        self.btn_export_plot = QPushButton("📸 导出视图为图片")
        self.btn_export_plot.clicked.connect(self.export_plot)

        self.btn_open_report = QPushButton("📄 打开纯文本计算书")
        self.btn_open_report.clicked.connect(self.open_text_report)

        ctrl_layout.addRow("变形放大倍数", self.scale_spin)
        ctrl_layout.addRow("反力箭头比例", self.reaction_scale_spin)
        ctrl_layout.addRow(self.chk_nodes)
        ctrl_layout.addRow(self.chk_elems)
        ctrl_layout.addRow(self.chk_deformed)
        ctrl_layout.addRow(self.chk_forces)
        ctrl_layout.addRow(self.chk_bmd)
        ctrl_layout.addRow(self.chk_reactions)
        ctrl_layout.addRow(self.btn_redraw)
        ctrl_layout.addRow(self.btn_export_plot)
        ctrl_layout.addRow(self.btn_open_report)

        self.viz_ctrl_dock.setWidget(ctrl_widget)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.viz_ctrl_dock)

    def _create_central_widget(self):
        central = QWidget()
        layout = QVBoxLayout()

        self.figure = Figure(figsize=(8, 6), dpi=110)
        self.figure.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.1)
        self.canvas = FigureCanvas(self.figure)
        self.nav_toolbar = NavigationToolbar(self.canvas, self)
        
        self.visualizer = FEMVisualizer(self.figure, self.canvas)
        self.visualizer.draw_empty()

        self.log_panel = QTextEdit()
        self.log_panel.setReadOnly(True)
        self.log_panel.setPlaceholderText("系统日志输出...")
        self.log_panel.setMaximumHeight(150)
        font = QFont("Consolas", 9)
        self.log_panel.setFont(font)

        layout.addWidget(self.nav_toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.log_panel)
        central.setLayout(layout)
        self.setCentralWidget(central)

    def _populate_tree(self):
        self.model_tree.clear()
        if not self.model:
            QTreeWidgetItem(self.model_tree, ["Nodes (0)"])
            QTreeWidgetItem(self.model_tree, ["Elements (0)"])
            QTreeWidgetItem(self.model_tree, ["Materials (0)"])
            QTreeWidgetItem(self.model_tree, ["Sections (0)"])
            QTreeWidgetItem(self.model_tree, ["Loads (0)"])
            QTreeWidgetItem(self.model_tree, ["Constraints (0)"])
            return

        nodes_item = QTreeWidgetItem(self.model_tree, [f"Nodes ({len(self.model.nodes)})"])
        for i, (ntype, x, y) in enumerate(self.model.nodes):
            t_str = "Frame" if ntype == 2 else "Truss"
            QTreeWidgetItem(nodes_item, [f"Node {i}: {t_str} ({x}, {y})"])

        elems_item = QTreeWidgetItem(self.model_tree, [f"Elements ({len(self.model.elements)})"])
        for i, (t, i0, i1, sec, mat) in enumerate(self.model.elements):
            QTreeWidgetItem(elems_item, [f"Elem {i}: N{i0}->N{i1} (Sec:{sec}, Mat:{mat})"])
            
        mats_item = QTreeWidgetItem(self.model_tree, [f"Materials ({len(self.model.materials)})"])
        for i, (e, mu, alpha) in enumerate(self.model.materials):
            QTreeWidgetItem(mats_item, [f"Mat {i}: E={e:g}, v={mu:g}"])

        sec_item = QTreeWidgetItem(self.model_tree, [f"Sections ({len(self.model.sections)})"])
        for i, (A, Iz, H) in enumerate(self.model.sections):
            QTreeWidgetItem(sec_item, [f"Sec {i}: A={A:g}, Iz={Iz:g}"])

        cons_item = QTreeWidgetItem(self.model_tree, [f"Constraints ({len(self.model.constraints)})"])
        for n, (cx, cy, cr) in self.model.constraints.items():
            fixed = []
            if cx < 0: fixed.append("X")
            if cy < 0: fixed.append("Y")
            if cr < 0: fixed.append("Rz")
            QTreeWidgetItem(cons_item, [f"Node {n}: Fixed({','.join(fixed)})"])

        loads_item = QTreeWidgetItem(self.model_tree, [f"Loads ({len(self.model.loads)})"])
        for i, (ltype, ldir, val, elem, node, pos, _, _) in enumerate(self.model.loads):
            dir_str = "X" if ldir == 0 else "Y" if ldir == 1 else "Rz"
            if ltype == 1:
                QTreeWidgetItem(loads_item, [f"Force: Node {node}, {dir_str}={val}"])
            else:
                QTreeWidgetItem(loads_item, [f"Load {i}: Type {ltype}, Val={val}"])

        self.model_tree.expandItem(nodes_item)
        self.model_tree.expandItem(elems_item)

    def _update_window_title(self):
        title = "FEM_2d GUI"
        if self.current_file:
            title += f" - {os.path.basename(self.current_file)}"
        if self.is_modified:
            title += " *"
        self.setWindowTitle(title)

    def _update_actions(self):
        has_model = self.model is not None
        self.save_action.setEnabled(has_model and self.is_modified)
        self.solve_action.setEnabled(has_model)
        self.export_plot_action.setEnabled(has_model and getattr(self.model, 'is_solved', False))

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
            show_bmd=self.chk_bmd.isChecked()
        )

    def _load_model_file(self, filepath):
        try:
            self.model = parse_input_file(filepath)
            self.current_file = filepath
            self.is_modified = False
            
            with open(filepath, 'r', encoding='utf-8') as f:
                self.text_editor.setPlainText(f.read())
            self.is_modified = False # Prevent setPlainText from triggering it
            
            self._populate_tree()
            self._update_window_title()
            self._update_actions()
            self._refresh_visualization()
            
            self.statusBar().showMessage(f"已加载: {os.path.basename(filepath)}", 3000)
            self._append_log(f"[INFO] 成功加载模型: {filepath}")
        except Exception as e:
            QMessageBox.critical(self, "加载失败", f"无法解析文件: {str(e)}")
            self._append_log(f"[ERROR] 读取失败: {str(e)}")

    @Slot()
    def _apply_editor_text(self):
        if not self.current_file:
            self.current_file = os.path.join(tempfile.gettempdir(), "temp_input.txt")
            
        try:
            with open(self.current_file, 'w', encoding='utf-8') as f:
                f.write(self.text_editor.toPlainText())
                
            self.model = parse_input_file(self.current_file)
            self.is_modified = False
            self._populate_tree()
            self._refresh_visualization()
            self._update_window_title()
            self._append_log("[INFO] 已从编辑器文字解析模型")
            self.statusBar().showMessage("模型解析成功", 3000)
        except Exception as e:
            QMessageBox.critical(self, "解析失败", f"文本格式有误:\\n{str(e)}")

    @Slot()
    def new_model(self):
        self.model = None
        self.current_file = None
        self.is_modified = False
        self.text_editor.clear()
        self.log_panel.clear()
        self.visualizer.draw_empty()
        self._populate_tree()
        self._update_window_title()
        self._update_actions()
        self.statusBar().showMessage("新建模型", 3000)

    @Slot()
    def open_default_example(self):
        example = os.path.join(self.app_dir, "example_continuous_beam.txt")
        if not os.path.exists(example):
            QMessageBox.warning(self, "文件丢失", f"找不到例子文件:\\n{example}")
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
            if not filename: return False
            self.current_file = filename
            
        try:
            with open(self.current_file, 'w', encoding='utf-8') as f:
                f.write(self.text_editor.toPlainText())
            self.is_modified = False
            self._update_window_title()
            self._update_actions()
            self.statusBar().showMessage("已保存", 3000)
            return True
        except Exception as e:
            QMessageBox.critical(self, "保存失败", f"无法写入文件:\\n{str(e)}")
            return False

    def _find_solver(self):
        """Find the fem_run executable."""
        import sys
        if getattr(sys, 'frozen', False):
            return os.path.join(self.app_dir, "bin", "fem_run.exe")
            
        candidates = [
            os.path.join(self.app_dir, "build", "bin", "Release", "fem_run.exe"),
            os.path.join(self.app_dir, "build", "bin", "Debug", "fem_run.exe"),
            os.path.join(self.app_dir, "build", "bin", "Release", "fem_run"),
        ]
        for c in candidates:
            if os.path.exists(c): return c
        return None

    @Slot()
    def solve_model(self):
        if not self.model:
            QMessageBox.warning(self, "无模型", "请先打开或创建模型")
            return
            
        if self.is_modified:
            self.save_model()

        exe_path = self._find_solver()
        if not exe_path:
            QMessageBox.critical(self, "找不到求解器", "未找到 fem_run.exe，请先编译项目(Release/Debug)")
            return

        os.makedirs(os.path.dirname(self.last_results_path), exist_ok=True)
        self.statusBar().showMessage("求解中...", 0)
        self._append_log(f"[INFO] 执行求解器: {exe_path}")

        try:
            solve_cmd = [exe_path, "--input", self.current_file, "--output", self.last_results_path, "--quiet"]
            
            # Using Popen for avoiding window flashes locally on Windows
            info = None
            if os.name == 'nt':
                info = subprocess.STARTUPINFO()
                info.dwFlags |= subprocess.STARTF_USESHOWWINDOW
                
            proc = subprocess.run(solve_cmd, cwd=self.work_dir, capture_output=True, text=True, startupinfo=info)
            
            if proc.stdout: self._append_log(proc.stdout.strip())
            if proc.stderr: self._append_log("[ERR] " + proc.stderr.strip())
            
            if proc.returncode != 0:
                QMessageBox.critical(self, "求解失败", "求解器执行失败，查看日志获取详情")
                self.statusBar().showMessage("求解失败", 5000)
                return

            self._append_log("[INFO] 解析结果文件...")
            self.model = parse_results_file(self.last_results_path, self.model)
            self._refresh_visualization()
            self._update_actions()
            self.statusBar().showMessage("求解完成！", 4000)

        except Exception as exc:
            self.statusBar().showMessage("发生错误", 5000)
            QMessageBox.critical(self, "错误", f"运行求解器时发生异常:\\n{exc}")

    @Slot()
    def export_plot(self):
        if not self.model or not self.model.is_solved:
            QMessageBox.information(self, "无结果", "请先求解模型再导出")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self, "导出图片", self.last_plot_path, "PNG 图片 (*.png);;全部文件 (*)"
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
            QMessageBox.information(self, "无计算书", "尚未在当前会话中成功求解，无法打开结果报告！")
            return
        try:
            if os.name == 'nt':
                os.startfile(self.last_results_path)
            else:
                import sys
                opener = "open" if sys.platform == "darwin" else "xdg-open"
                subprocess.run([opener, self.last_results_path])
        except Exception as e:
            QMessageBox.critical(self, "打开失败", f"无法用系统默认应用打开文件:\n{e}")

    @Slot()
    def show_about(self):
        QMessageBox.about(
            self, "关于 FEM_2d",
            "<h3>FEM_2d</h3><p>版本 2.0 (已重构)</p><p>二维杆系与框架有限元分析程序GUI</p>"
        )

    def closeEvent(self, event):
        if self.is_modified:
            reply = QMessageBox.question(
                self, "未保存", "有未保存的修改，是否在退出前保存？",
                QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel
            )
            if reply == QMessageBox.StandardButton.Save:
                if not self.save_model(): event.ignore(); return
            elif reply == QMessageBox.StandardButton.Cancel:
                event.ignore(); return
        event.accept()
