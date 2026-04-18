"""
Main Window for FEM_2d GUI Application

This module implements the main application window with menu bar, toolbar,
status bar, and docked panels for model tree and properties.
"""

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QMenuBar, QMenu, QToolBar, QStatusBar, QDockWidget,
    QTreeWidget, QTreeWidgetItem, QTextEdit, QLabel,
    QFileDialog, QMessageBox, QSplitter
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QAction, QIcon, QKeySequence

import os


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
        self.solve_action.setStatusTip("Run FEM analysis")
        self.solve_action.triggered.connect(self.solve_model)

        # View actions
        self.show_deformed_action = QAction("Show &Deformed", self)
        self.show_deformed_action.setShortcut(Qt.Key_F6)
        self.show_deformed_action.setCheckable(True)

        self.show_forces_action = QAction("Show &Forces", self)
        self.show_forces_action.setShortcut(Qt.Key_F7)
        self.show_forces_action.setCheckable(True)

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

        # View menu
        view_menu = menubar.addMenu("&View")
        view_menu.addAction(self.show_deformed_action)
        view_menu.addAction(self.show_forces_action)

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
        main_toolbar.addAction(self.save_action)
        main_toolbar.addSeparator()
        main_toolbar.addAction(self.solve_action)

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

    def _create_central_widget(self):
        """Create central widget with visualization."""
        central = QWidget()
        layout = QVBoxLayout()

        # Placeholder for visualization
        self.viz_placeholder = QLabel("Visualization Area\n\n(PyVista view will be integrated here)")
        self.viz_placeholder.setAlignment(Qt.AlignCenter)
        self.viz_placeholder.setStyleSheet("""
            QLabel {
                background-color: #f0f0f0;
                border: 2px dashed #cccccc;
                font-size: 16px;
                color: #999999;
            }
        """)

        layout.addWidget(self.viz_placeholder)
        central.setLayout(layout)
        self.setCentralWidget(central)

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
        has_model = self.model is not None
        self.save_action.setEnabled(has_model and self.is_modified)
        self.save_as_action.setEnabled(has_model)
        self.solve_action.setEnabled(has_model)

    # Slots for actions
    @Slot()
    def new_model(self):
        """Create a new model."""
        # TODO: Check for unsaved changes
        self.model = None  # Will be FEM Model
        self.current_file = None
        self.is_modified = False
        self._populate_tree()
        self._update_window_title()
        self._update_actions()
        self.statusBar().showMessage("New model created", 3000)

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
            # TODO: Load model
            self.current_file = filename
            self.is_modified = False
            self._update_window_title()
            self._update_actions()
            self.statusBar().showMessage(f"Loaded: {os.path.basename(filename)}", 3000)

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
        """Solve the FEM problem."""
        if not self.model:
            QMessageBox.warning(self, "No Model", "Please create or load a model first.")
            return

        self.statusBar().showMessage("Solving...", 0)
        # TODO: Run solver
        self.statusBar().showMessage("Solution complete", 3000)

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
