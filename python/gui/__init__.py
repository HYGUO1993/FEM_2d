"""
GUI package for FEM_2d application

This package contains all GUI-related modules including the main window,
visualization widgets, property editors, and dialogs.
"""

from .main_window import MainWindow
from .theme import apply_theme, setup_matplotlib_dark_theme
from .fem_parser import parse_input_file, parse_results_file, FEMModelData
from .visualization import FEMVisualizer

__all__ = [
    "MainWindow", 
    "apply_theme", 
    "setup_matplotlib_dark_theme",
    "parse_input_file",
    "parse_results_file",
    "FEMModelData",
    "FEMVisualizer"
]
