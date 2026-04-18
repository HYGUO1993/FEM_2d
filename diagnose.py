#!/usr/bin/env python3
"""
FEM_2d Environment Diagnostic Tool

Use this script to diagnose environment and dependency issues.
"""

import sys
import os
import platform

def check_python():
    """Check Python version and executable."""
    print("=" * 60)
    print("PYTHON ENVIRONMENT")
    print("=" * 60)
    print(f"Python version:   {sys.version}")
    print(f"Executable:       {sys.executable}")
    print(f"Platform:         {platform.platform()}")
    print(f"Architecture:     {platform.architecture()}")
    print()

def check_packages():
    """Check installed packages."""
    print("=" * 60)
    print("PACKAGE CHECK")
    print("=" * 60)
    
    packages = {
        'PySide6': 'GUI Framework',
        'matplotlib': 'Plotting Library',
        'numpy': 'Numerical Computing',
    }
    
    for pkg_name, description in packages.items():
        try:
            mod = __import__(pkg_name)
            version = getattr(mod, '__version__', 'unknown')
            print("[OK] {:<15} {:<20} ({})".format(pkg_name, version, description))
        except ImportError as e:
            print("[NO] {:<15} NOT INSTALLED       ({})".format(pkg_name, description))
    
    print()

def check_fem_backend():
    """Check FEM C++ backend."""
    print("=" * 60)
    print("FEM SOLVER")
    print("=" * 60)
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    fem_exe_windows = os.path.join(project_root, "build", "bin", "Release", "fem_run.exe")
    fem_exe_unix = os.path.join(project_root, "build", "bin", "Release", "fem_run")
    
    if os.path.exists(fem_exe_windows):
        print("[OK] Solver (Windows): {}".format(fem_exe_windows))
    elif os.path.exists(fem_exe_unix):
        print("[OK] Solver (Unix):    {}".format(fem_exe_unix))
    else:
        print("[NO] Solver not found - Build with:")
        print("    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release")
        print("    cmake --build build --config Release")
    
    print()

def check_paths():
    """Check Python paths and project structure."""
    print("=" * 60)
    print("PROJECT STRUCTURE")
    print("=" * 60)
    
    project_root = os.path.dirname(os.path.abspath(__file__))
    print(f"Project root:     {project_root}")
    
    required_dirs = {
        'python': 'Python GUI code',
        'python/gui': 'GUI modules',
        'build': 'Build output',
    }
    
    for dirname, desc in required_dirs.items():
        path = os.path.join(project_root, dirname)
        exists = "[OK]" if os.path.exists(path) else "[NO]"
        print("{} {:<20} ({})".format(exists, dirname, desc))
    
    print()

def test_imports():
    """Test if critical imports work."""
    print("=" * 60)
    print("IMPORT TEST")
    print("=" * 60)
    
    issues = []
    
    # Test PySide6 import with details
    pyside6_ok = False
    try:
        from PySide6.QtWidgets import QApplication
        print("[OK] PySide6.QtWidgets imported successfully")
        pyside6_ok = True
    except ModuleNotFoundError as e:
        print("[NO] PySide6 not installed: {}".format(e))
        issues.append(('pyside6', 'notinstalled', str(e)))
    except (ImportError, OSError) as e:
        error_msg = str(e)
        print("[NO] Error loading PySide6: {}".format(error_msg))
        
        # Check for DLL loading errors on Windows
        if "DLL load failed" in error_msg or "procedure entry point" in error_msg or "找不到指定的程序" in error_msg:
            issues.append(('pyside6', 'dll', error_msg))
        else:
            issues.append(('pyside6', 'runtime', error_msg))
    except Exception as e:
        error_msg = "{}: {}".format(type(e).__name__, e)
        print("[NO] Unexpected error: {}".format(error_msg))
        issues.append(('pyside6', 'unexpected', error_msg))
    
    try:
        import matplotlib.pyplot
        print("[OK] matplotlib.pyplot imported successfully")
    except ImportError as e:
        print("[NO] Failed to import matplotlib: {}".format(e))
        issues.append(('matplotlib', 'import', str(e)))
    except Exception as e:
        print("[NO] Error loading matplotlib: {}: {}".format(type(e).__name__, e))
        issues.append(('matplotlib', 'runtime', str(e)))
    
    try:
        import numpy
        print("[OK] numpy imported successfully")
    except ImportError as e:
        print("[NO] Failed to import numpy: {}".format(e))
        issues.append(('numpy', 'import', str(e)))
    except Exception as e:
        print("[NO] Error loading numpy: {}: {}".format(type(e).__name__, e))
        issues.append(('numpy', 'runtime', str(e)))
    
    print()
    return issues

def suggest_fixes(import_issues):
    """Suggest fixes based on detected issues."""
    print("=" * 60)
    print("TROUBLESHOOTING")
    print("=" * 60)
    
    if not import_issues:
        print("[OK] No obvious issues detected!")
        print()
        print("If GUI still won't launch, try:")
        print("  1. Run: python quick_start_gui.py")
        print("  2. Or: python python/fem_gui.py")
    else:
        print("[!!] Issues detected:")
        print()
        
        # Check for DLL errors (Windows specific)
        has_dll_error = any(issue[1] == 'dll' for issue in import_issues)
        
        if has_dll_error:
            print("[!!] WINDOWS DLL LOADING ERROR DETECTED")
            print()
            print("This means PySide6 package is installed, but Windows is missing")
            print("required system libraries (Visual C++ runtime, Qt libraries, etc.)")
            print()
            print("Solutions (try in order):")
            print("  1. Install Visual C++ Redistributable:")
            print("     https://support.microsoft.com/en-us/help/2977003")
            print("     -> Download & install 'vc_redist.x64.exe' (or x86 if 32-bit Python)")
            print()
            print("  2. Reinstall PySide6 after installing VC++ runtime:")
            print("     python -m pip cache purge")
            print("     python -m pip install --force-reinstall PySide6")
            print()
            print("  3. Or use conda (better system library handling):")
            print("     conda install -c conda-forge PySide6")
            print()
            print("[!] Or run auto-fix: fix_dll_error.bat (or 修复DLL错误.bat)")
            print()
            print("[!] See WINDOWS_DLL_ERROR.md for complete solutions")
        
        # Check for missing packages
        missing_packages = [issue[0] for issue in import_issues if issue[1] == 'notinstalled']
        if missing_packages:
            print()
            print("Missing packages: {}".format(", ".join(missing_packages)))
            print("Install with: python -m pip install {}".format(" ".join(missing_packages)))
    
    print()

def _is_installed(pkg_name):
    """Check if a package is installed."""
    try:
        __import__(pkg_name)
        return True
    except ImportError:
        return False

def main():
    """Run all diagnostics."""
    print()
    print("╔" + "=" * 58 + "╗")
    print("║" + " FEM_2d Environment Diagnostic Tool".center(58) + "║")
    print("╚" + "=" * 58 + "╝")
    print()
    
    check_python()
    check_packages()
    check_fem_backend()
    check_paths()
    import_issues = test_imports()
    suggest_fixes(import_issues)
    
    print("=" * 60)
    print("For more help, see:")
    print("  - QUICKSTART.md")
    print("  - GUI_USAGE.md")
    print("  - INSTALL_DEPENDENCIES.md")
    print("  - WINDOWS_DLL_ERROR.md (if you see DLL errors)")
    print("  - 修复依赖问题.md")
    print("=" * 60)

if __name__ == "__main__":
    main()
