@echo off
REM FEM_2d GUI Launcher (Windows Batch)
REM 
REM This script launches the FEM_2d GUI application.
REM It checks for required dependencies and handles environment setup.
REM
REM Usage: 
REM   - Double-click this file to launch the GUI
REM   - Or run: cmd /c scripts\run_gui.bat
REM
REM Requirements:
REM   - Python 3.8+ with conda or system Python
REM   - PySide6, matplotlib, numpy packages

setlocal enabledelayedexpansion

REM Get the project root directory
cd /d "%~dp0.."
set PROJECT_ROOT=%cd%

REM Check if Python is available
where python >nul 2>&1
if errorlevel 1 (
    echo Error: Python not found in PATH
    echo Please install Python 3.8+ or activate a conda environment with Python
    echo.
    pause
    exit /b 1
)

echo ============================================================
echo FEM_2d GUI Application Launcher
echo ============================================================
echo.

REM Check Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo [INFO] Python version: %PYTHON_VERSION%
echo [INFO] Project root: %PROJECT_ROOT%
echo.

REM Check dependencies
echo [INFO] Checking dependencies...
set MISSING_PACKAGES=

python -c "import PySide6" 2>nul
if errorlevel 1 set MISSING_PACKAGES=%MISSING_PACKAGES% PySide6

python -c "import matplotlib" 2>nul
if errorlevel 1 set MISSING_PACKAGES=%MISSING_PACKAGES% matplotlib

python -c "import numpy" 2>nul
if errorlevel 1 set MISSING_PACKAGES=%MISSING_PACKAGES% numpy

if not "!MISSING_PACKAGES!"=="" (
    echo.
    echo [WARNING] Missing packages:!MISSING_PACKAGES!
    echo.
    echo [INFO] Installing missing packages...
    pip install -q --no-warn-script-location PySide6 matplotlib numpy
    if errorlevel 1 (
        echo [ERROR] Failed to install dependencies
        echo Please install manually with:
        echo   pip install PySide6 matplotlib numpy
        echo.
        pause
        exit /b 1
    )
)

echo [INFO] All dependencies available
echo.

REM Check if FEM release executable is built
set FEM_EXE=%PROJECT_ROOT%\build\bin\Release\fem_run.exe
if not exist "%FEM_EXE%" (
    echo [WARNING] Release executable not found: %FEM_EXE%
    echo [INFO] The solver will not work until you build the Release version:
    echo   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    echo   cmake --build build --config Release
    echo.
)

echo [INFO] Launching GUI...
echo ============================================================
echo.

REM Launch the GUI
cd /d "%PROJECT_ROOT%"
python python\fem_gui.py

REM Cleanup on exit
echo.
echo [INFO] GUI closed
exit /b 0
