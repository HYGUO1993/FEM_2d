from cx_Freeze import setup, Executable
import sys
import os

sys.path.append(os.path.abspath("python"))

build_exe_options = {
    "packages": ["gui", "os", "sys", "PyQt6", "matplotlib", "numpy"],
    "include_files": [
        ("build/bin/Release/fem_run.exe", "bin/fem_run.exe"),
        ("example_continuous_beam.txt", "example_continuous_beam.txt"),
        ("example_multi_story_frame.txt", "example_multi_story_frame.txt"),
        ("example_truss_bridge.txt", "example_truss_bridge.txt")
    ],
    "excludes": ["tkinter"],
}

base = "Win32GUI" if sys.platform == "win32" else None

setup(
    name="FEM_2d_Studio",
    version="1.0",
    description="FEM_2d GUI Application",
    options={"build_exe": build_exe_options},
    executables=[Executable("python/fem_gui.py", base=base, target_name="FEM_2d_Studio.exe")]
)
