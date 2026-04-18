"""
FEM_2d - 2D Finite Element Analysis Library
=============================================

A modern finite element analysis library for 2D truss and frame structures,
designed for engineering education and practice.

Basic Usage
-----------

Using the C API (via ctypes):
    >>> from fem2d import FEMModel
    >>> model = FEMModel()
    >>> model.add_material(E=2.1e11, mu=0.3, alpha=1.2e-5)
    >>> model.add_section(A=0.01, Iz=8.333e-6, H=0.1)
    >>> # Add nodes, elements, constraints, loads...
    >>> model.solve()
    >>> displacements = model.get_displacements()

Using pybind11 bindings (if available):
    >>> import pyfem
    >>> model = pyfem.Model()
    >>> model.add_material(2.1e11, 0.3, 1.2e-5)
    >>> # ... similar to above

Modules
-------
- model: High-level Python interface to FEM solver
- visualization: Plotting and visualization utilities
- examples: Example problems and tutorials
"""

__version__ = "1.0.0"
__author__ = "FEM_2d Development Team"

# Try to import pybind11 version first, fall back to ctypes
try:
    from . import pyfem
    _BACKEND = "pybind11"
    Model = pyfem.Model
except ImportError:
    try:
        from .model_ctypes import FEMModel as Model
        _BACKEND = "ctypes"
    except ImportError:
        # If neither is available, create a dummy class
        class Model:
            def __init__(self):
                raise ImportError(
                    "FEM backend not available. Please build the C++ library first."
                )
        _BACKEND = "none"

# Import visualization tools
try:
    from .visualization import visualize_results, plot_structure
except ImportError:
    pass

__all__ = [
    "Model",
    "visualize_results",
    "plot_structure",
    "__version__",
]
