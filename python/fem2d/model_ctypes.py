"""
FEM Model interface using ctypes (C API)

This module provides a Python interface to the FEM C API using ctypes.
It's used when pybind11 bindings are not available.
"""

import ctypes
import os
import sys
import numpy as np
from typing import Optional, Tuple


# Constants
TRUSS = 1
FRAME = 2
TRUSS_NODE = 1
FRAME_NODE = 2
FORCE_ON_NODE = 1
DIRECT_X = 0
DIRECT_Y = 1
DIRECT_R = 2


def _find_library():
    """Find the FEM shared library."""
    possible_names = [
        "libfem_core.so",      # Linux
        "libfem_core.dylib",   # macOS
        "fem_core.dll",        # Windows
        "libfem_core.dll",     # Windows alternative
    ]

    # Search paths
    search_paths = [
        os.path.join(os.path.dirname(__file__), "..", "..", "build", "lib"),
        os.path.join(os.path.dirname(__file__), "..", "..", "build", "bin"),
        "/usr/local/lib",
        "/usr/lib",
    ]

    for path in search_paths:
        for name in possible_names:
            full_path = os.path.join(path, name)
            if os.path.exists(full_path):
                return full_path

    raise RuntimeError(
        f"Could not find FEM shared library. Searched in: {search_paths}"
    )


# Load library
try:
    _lib_path = _find_library()
    _lib = ctypes.CDLL(_lib_path)
except Exception as e:
    _lib = None
    print(f"Warning: Could not load FEM library: {e}", file=sys.stderr)


def _setup_api():
    """Setup function signatures for C API."""
    if _lib is None:
        return

    # FEM_CreateModel
    _lib.FEM_CreateModel.argtypes = []
    _lib.FEM_CreateModel.restype = ctypes.c_void_p

    # FEM_DestroyModel
    _lib.FEM_DestroyModel.argtypes = [ctypes.c_void_p]
    _lib.FEM_DestroyModel.restype = None

    # FEM_AddNode
    _lib.FEM_AddNode.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                   ctypes.c_double, ctypes.c_double]
    _lib.FEM_AddNode.restype = ctypes.c_int

    # FEM_AddElement
    _lib.FEM_AddElement.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                      ctypes.c_int, ctypes.c_int,
                                      ctypes.c_int, ctypes.c_int]
    _lib.FEM_AddElement.restype = ctypes.c_int

    # FEM_AddMaterial
    _lib.FEM_AddMaterial.argtypes = [ctypes.c_void_p, ctypes.c_double,
                                       ctypes.c_double, ctypes.c_double]
    _lib.FEM_AddMaterial.restype = ctypes.c_int

    # FEM_AddSection
    _lib.FEM_AddSection.argtypes = [ctypes.c_void_p, ctypes.c_double,
                                      ctypes.c_double, ctypes.c_double]
    _lib.FEM_AddSection.restype = ctypes.c_int

    # FEM_AddLoad
    _lib.FEM_AddLoad.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                                   ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                   ctypes.c_double, ctypes.c_double, ctypes.c_double]
    _lib.FEM_AddLoad.restype = ctypes.c_int

    # FEM_AddConstraint
    _lib.FEM_AddConstraint.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                         ctypes.c_int, ctypes.c_int, ctypes.c_int]
    _lib.FEM_AddConstraint.restype = ctypes.c_int

    # FEM_LoadFromFile
    _lib.FEM_LoadFromFile.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    _lib.FEM_LoadFromFile.restype = ctypes.c_int

    # FEM_Solve
    _lib.FEM_Solve.argtypes = [ctypes.c_void_p]
    _lib.FEM_Solve.restype = ctypes.c_int

    # FEM_ExportResults
    _lib.FEM_ExportResults.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    _lib.FEM_ExportResults.restype = ctypes.c_int

    # FEM_GetDisplacements
    _lib.FEM_GetDisplacements.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
    _lib.FEM_GetDisplacements.restype = ctypes.POINTER(ctypes.c_double)

    # FEM_GetNodeCount
    _lib.FEM_GetNodeCount.argtypes = [ctypes.c_void_p]
    _lib.FEM_GetNodeCount.restype = ctypes.c_int

    # FEM_GetLastError
    _lib.FEM_GetLastError.argtypes = []
    _lib.FEM_GetLastError.restype = ctypes.c_char_p


_setup_api()


class FEMModel:
    """Python wrapper for FEM C API."""

    def __init__(self):
        """Create a new FEM model."""
        if _lib is None:
            raise RuntimeError("FEM library not loaded")
        self._handle = _lib.FEM_CreateModel()
        if not self._handle:
            raise RuntimeError("Failed to create FEM model")

    def __del__(self):
        """Destroy the model."""
        if hasattr(self, '_handle') and self._handle:
            _lib.FEM_DestroyModel(self._handle)

    def add_node(self, node_type: int, x: float, y: float) -> int:
        """Add a node to the model."""
        return _lib.FEM_AddNode(self._handle, node_type, x, y)

    def add_element(self, elem_type: int, node1: int, node2: int,
                   section: int, material: int) -> int:
        """Add an element to the model."""
        return _lib.FEM_AddElement(self._handle, elem_type, node1, node2,
                                   section, material)

    def add_material(self, E: float, mu: float, alpha: float) -> int:
        """Add a material to the model."""
        return _lib.FEM_AddMaterial(self._handle, E, mu, alpha)

    def add_section(self, A: float, Iz: float, H: float) -> int:
        """Add a section to the model."""
        return _lib.FEM_AddSection(self._handle, A, Iz, H)

    def add_load(self, load_type: int, direction: int, value: float,
                elem: int = -1, node: int = -1, position: float = 0.0,
                T0: float = 0.0, T1: float = 0.0) -> int:
        """Add a load to the model."""
        return _lib.FEM_AddLoad(self._handle, load_type, direction, value,
                               elem, node, position, T0, T1)

    def add_constraint(self, node: int, dof_x: int, dof_y: int, dof_r: int) -> int:
        """Add a constraint to the model."""
        return _lib.FEM_AddConstraint(self._handle, node, dof_x, dof_y, dof_r)

    def load_from_file(self, filename: str) -> bool:
        """Load model from file."""
        result = _lib.FEM_LoadFromFile(self._handle, filename.encode('utf-8'))
        return result == 0

    def solve(self) -> bool:
        """Solve the FEM problem."""
        result = _lib.FEM_Solve(self._handle)
        return result == 0

    def export_results(self, filename: str) -> bool:
        """Export results to file."""
        result = _lib.FEM_ExportResults(self._handle, filename.encode('utf-8'))
        return result == 0

    def get_displacements(self) -> Optional[np.ndarray]:
        """Get displacements as numpy array."""
        count = ctypes.c_int()
        ptr = _lib.FEM_GetDisplacements(self._handle, ctypes.byref(count))
        if not ptr or count.value == 0:
            return None
        return np.ctypeslib.as_array(ptr, shape=(count.value,))

    @property
    def node_count(self) -> int:
        """Get number of nodes."""
        return _lib.FEM_GetNodeCount(self._handle)

    def get_last_error(self) -> str:
        """Get last error message."""
        err = _lib.FEM_GetLastError()
        return err.decode('utf-8') if err else ""
