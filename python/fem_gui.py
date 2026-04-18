#!/usr/bin/env python3
"""
FEM_2d GUI Application

A modern graphical user interface for 2D finite element analysis,
designed for engineering education and practice.

Requirements:
    - PySide6 (Qt for Python)
    - Matplotlib (2D plotting)
    - NumPy

Installation:
    pip install PySide6 matplotlib numpy
"""

import sys
import os
import traceback

# Add parent directory to path to import fem2d
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

def try_import_pyside6():
    """Try to import PySide6 with detailed error handling."""
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
        from PySide6.QtCore import Qt
        return True, None
    except ImportError as e:
        return False, str(e)
    except Exception as e:
        # Other errors (like OpenGL, library loading issues)
        return False, f"Runtime error: {str(e)}"

pyside6_ok, pyside6_error = try_import_pyside6()

if not pyside6_ok:
    print("Error: Failed to load PySide6")
    print(f"Details: {pyside6_error}")
    print()
    print("Solutions:")
    print("  1. Install PySide6: pip install --user PySide6")
    print("  2. Or use conda: conda install -c conda-forge PySide6")
    print()
    if "PySide6" in pyside6_error or "No module named" in pyside6_error:
        print("PySide6 package not found. Please install it.")
    else:
        print("PySide6 is installed but failed to load (may need system libraries).")
    sys.exit(1)

from PySide6.QtWidgets import QApplication, QMessageBox
from PySide6.QtCore import Qt
from gui.main_window import MainWindow


def check_dependencies():
    """Check if all required dependencies are available."""
    missing = []

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
