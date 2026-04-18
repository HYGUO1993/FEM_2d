@echo off
REM ============================================================
REM FEM_2d GUI Launcher (Simplified)
REM 在项目根目录双击即可启动 GUI
REM ============================================================

setlocal enabledelayedexpansion

REM Get the project root directory
set PROJECT_ROOT=%~dp0

REM Change to project root
cd /d "%PROJECT_ROOT%"

REM Check if Python is available
where python >nul 2>&1
if errorlevel 1 (
    echo.
    echo [ERROR] Python not found in PATH
    echo.
    echo 请安装 Python 3.8+ 或激活 conda 环境后重试
    echo.
    pause
    exit /b 1
)

REM Install dependencies if needed
echo [INFO] Checking dependencies...
python -c "import PySide6; import matplotlib; import numpy" 2>nul
if errorlevel 1 (
    echo.
    echo [INFO] Installing required packages...
    echo.
    python -m pip install --upgrade pip -q 2>nul
    python -m pip install PySide6 matplotlib numpy
    if errorlevel 1 (
        echo.
        echo [ERROR] Failed to install dependencies automatically
        echo.
        echo Please try to install manually by running:
        echo   python -m pip install --user PySide6 matplotlib numpy
        echo.
        echo Or using conda:
        echo   conda install -c conda-forge PySide6 matplotlib numpy
        echo.
        pause
        exit /b 1
    )
    echo [INFO] Dependencies installed successfully
)

REM Launch GUI
echo [INFO] Starting FEM_2d GUI...
python python\fem_gui.py

REM Cleanup
exit /b 0
