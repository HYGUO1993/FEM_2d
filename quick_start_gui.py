#!/usr/bin/env python3
"""
Quick Start: Launch FEM_2d GUI

This script provides a simple way to start the FEM_2d GUI application
with proper error handling and user feedback.

Usage:
    python quick_start_gui.py
"""

import sys
import os
import subprocess
import platform

def main():
    """Main launcher function."""
    print("=" * 60)
    print("FEM_2d GUI - Quick Start Launcher")
    print("=" * 60)
    print()
    
    # Get project root
    project_root = os.path.dirname(os.path.abspath(__file__))
    gui_script = os.path.join(project_root, "python", "fem_gui.py")
    
    # Check if GUI script exists
    if not os.path.exists(gui_script):
        print(f"[ERROR] GUI script not found: {gui_script}")
        return 1
    
    print(f"[INFO] Project root: {project_root}")
    print(f"[INFO] Python version: {sys.version.split()[0]}")
    print()
    
    # Check Python version
    if sys.version_info < (3, 8):
        print(f"[ERROR] Python 3.8+ required, got {sys.version}")
        return 1
    
    # Check dependencies
    print("[INFO] Checking dependencies...")
    missing = []
    
    try:
        import PySide6
        print("  ✓ PySide6")
    except ImportError:
        print("  ✗ PySide6 (missing)")
        missing.append("PySide6")
    
    try:
        import matplotlib
        print("  ✓ matplotlib")
    except ImportError:
        print("  ✗ matplotlib (missing)")
        missing.append("matplotlib")
    
    try:
        import numpy
        print("  ✓ numpy")
    except ImportError:
        print("  ✗ numpy (missing)")
        missing.append("numpy")
    
    if missing:
        print()
        print("[WARNING] Missing packages:", ", ".join(missing))
        print("[INFO] Installing missing packages...")
        print()
        
        # Try with --user flag if regular install fails
        install_methods = [
            [sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
            [sys.executable, "-m", "pip", "install"] + missing,
            [sys.executable, "-m", "pip", "install", "--user"] + missing,
        ]
        
        success = False
        for method in install_methods:
            try:
                subprocess.check_call(method, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                print("[INFO] Installation successful")
                success = True
                break
            except subprocess.CalledProcessError:
                continue
        
        if not success:
            print("[ERROR] Failed to install dependencies automatically")
            print()
            print("Please install manually using one of these methods:")
            print()
            print("Method 1 (pip with --user flag):")
            print(f"  python -m pip install --user {' '.join(missing)}")
            print()
            print("Method 2 (using conda):")
            print(f"  conda install -c conda-forge {' '.join(missing)}")
            print()
            print("Method 3 (pip with sudo on Linux/macOS):")
            print(f"  sudo pip install {' '.join(missing)}")
            print()
            return 1
    
    print()
    
    # Check FEM solver
    fem_exe = os.path.join(project_root, "build", "bin", "Release", "fem_run.exe")
    if platform.system() != "Windows":
        fem_exe = os.path.join(project_root, "build", "bin", "Release", "fem_run")
    
    if os.path.exists(fem_exe):
        print(f"[INFO] FEM solver found: {fem_exe}")
    else:
        print("[WARNING] FEM solver executable not found")
        print("[INFO] To enable solver, build the Release version:")
        print("    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release")
        print("    cmake --build build --config Release")
    
    print()
    print("[INFO] Launching GUI...")
    print("=" * 60)
    print()
    
    # Launch GUI
    try:
        subprocess.run([sys.executable, gui_script], cwd=project_root, check=False)
    except Exception as e:
        print(f"[ERROR] Failed to launch GUI: {e}")
        return 1
    
    print()
    print("[INFO] GUI closed")
    return 0

if __name__ == "__main__":
    sys.exit(main())
