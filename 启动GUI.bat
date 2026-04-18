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

REM Activate conda environment
echo [INFO] Activating conda environment 'fem_gui_env'...
call conda activate fem_gui_env 2>nul
if errorlevel 1 (
    echo [INFO] Environment 'fem_gui_env' not found.
    echo [INFO] Creating 'fem_gui_env' with required dependencies...
    call conda create -y -n fem_gui_env python=3.10 PySide6 matplotlib numpy -c conda-forge
    call conda activate fem_gui_env 2>nul
    if errorlevel 1 (
        echo.
        echo [ERROR] Failed to activate conda environment.
        echo 请确保已正确安装 Anaconda/Miniconda 并将其添加到环境变量中。
        echo.
        pause
        exit /b 1
    )
)

REM Launch GUI
echo [INFO] Starting FEM_2d GUI...
python python\fem_gui.py

REM Cleanup
exit /b 0
