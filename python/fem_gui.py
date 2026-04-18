#!/usr/bin/env python3
"""
FEM_2d GUI Application

A modern graphical user interface for 2D finite element analysis,
designed for engineering education and practice.

Requirements:
    - PySide6 (Qt for Python)
    - PyVista (3D visualization)
    - Matplotlib (2D plotting)
    - NumPy

Installation:
    pip install PySide6 pyvista matplotlib numpy
"""

import sys
import os

# Add parent directory to path to import fem2d
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    from PySide6.QtWidgets import QApplication, QMessageBox
    from PySide6.QtCore import Qt
except ImportError:
    print("Error: PySide6 is not installed.")
    print("Please install it with: pip install PySide6")
    sys.exit(1)

from gui.main_window import MainWindow


def check_dependencies():
    """Check if all required dependencies are available."""
    missing = []

    # Check for PyVista
    try:
        import pyvista
    except ImportError:
        missing.append("pyvista")

    # Check for Matplotlib
    try:
        import matplotlib
    except ImportError:
        missing.append("matplotlib")

    # Check for NumPy
    try:
        import numpy
    except ImportError:
        missing.append("numpy")

    # Check for FEM backend
    try:
        from fem2d import Model
    except ImportError:
        print("Warning: FEM backend not available.")
        print("Some features will be limited until the C++ library is built.")

    if missing:
        msg = "The following required packages are missing:\n\n"
        msg += "\n".join(f"  - {pkg}" for pkg in missing)
        msg += "\n\nPlease install them with:\n"
        msg += f"pip install {' '.join(missing)}"
        return msg

    return None


def main():
    """Main application entry point."""
    # Create application
    app = QApplication(sys.argv)
    app.setApplicationName("FEM_2d")
    app.setOrganizationName("FEM_2d Development Team")
    app.setApplicationVersion("1.0.0")

    # Set application style
    app.setStyle("Fusion")

    # Check dependencies
    dep_error = check_dependencies()
    if dep_error:
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setWindowTitle("Missing Dependencies")
        msg_box.setText("Some required packages are missing.")
        msg_box.setDetailedText(dep_error)
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.exec()
        # Continue anyway, some features may work

    # Create and show main window
    window = MainWindow()
    window.show()

    # Run application
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
