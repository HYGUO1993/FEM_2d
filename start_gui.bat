@echo off
REM ============================================================
REM FEM_2d GUI Launcher - Double-click to run!
REM ============================================================

setlocal enabledelayedexpansion

set PROJECT_ROOT=%~dp0
cd /d "%PROJECT_ROOT%"

echo.
echo ============================================================
echo FEM_2d GUI - Starting...
echo ============================================================
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found!
    echo Please install Python 3.8+ or activate conda
    echo.
    pause
    exit /b 1
)

echo [1/2] Checking dependencies...
python -c "import PySide6, matplotlib, numpy" 2>nul
if errorlevel 1 (
    echo [2/2] Installing required packages...
    python -m pip install --upgrade pip -q 2>nul
    python -m pip install PySide6 matplotlib numpy
    if errorlevel 1 (
        echo.
        echo [ERROR] Failed to install dependencies
        echo.
        echo Try: python -m pip install --user PySide6 matplotlib numpy
        echo.
        pause
        exit /b 1
    )
)

echo [2/2] Launching GUI...
python python\fem_gui.py
exit /b 0
