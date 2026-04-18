@echo off
REM ============================================================
REM Install FEM_2d GUI Dependencies
REM 用于修复依赖安装问题
REM ============================================================

setlocal enabledelayedexpansion

echo.
echo ============================================================
echo FEM_2d GUI Dependencies Installer
echo ============================================================
echo.

REM Check Python
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found!
    echo Please install Python 3.8+ first
    pause
    exit /b 1
)

echo [1/3] Upgrading pip...
python -m pip install --upgrade pip --quiet
if errorlevel 1 echo [WARNING] pip upgrade failed, continuing anyway...

echo [2/3] Installing required packages...
echo.

REM Try multiple installation methods
python -m pip install PySide6 matplotlib numpy
if errorlevel 1 (
    echo.
    echo [WARNING] Standard installation failed, trying --user flag...
    python -m pip install --user PySide6 matplotlib numpy
    if errorlevel 1 (
        echo.
        echo [ERROR] Installation failed with both methods
        echo.
        echo Try using conda instead:
        echo   conda install -c conda-forge PySide6 matplotlib numpy
        echo.
        pause
        exit /b 1
    )
)

echo.
echo [3/3] Verifying installation...
python -c "import PySide6; import matplotlib; import numpy; print('[SUCCESS] All packages installed!')"
if errorlevel 1 (
    echo [ERROR] Verification failed
    pause
    exit /b 1
)

echo.
echo ============================================================
echo [SUCCESS] Dependencies installed successfully!
echo ============================================================
echo.
echo You can now run:
echo   - Double-click: 启动GUI.bat or start_gui.bat
echo   - Command: python quick_start_gui.py
echo   - VS Code: Press Ctrl+Shift+P and search "Launch FEM_2d GUI"
echo.
pause
exit /b 0
