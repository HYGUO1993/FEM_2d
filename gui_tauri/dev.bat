@echo off
REM ============================================================
REM FemLab Studio (Tauri) - 开发模式一键启动
REM 直接运行编译好的 exe (frontendDist 已嵌入, 无需外部 dev server)
REM 需要: 已编译的 femlab-studio.exe + femcli.exe
REM 首次: 先运行 build_gui.bat 编译 Tauri 后端
REM ============================================================
setlocal
cd /d "%~dp0"

echo ============================================================
echo  FemLab Studio (Tauri) - 开发模式
echo ============================================================

REM 检查 femcli
if not exist "..\build\bin\femcli.exe" (
    echo [错误] 未找到 ..\build\bin\femcli.exe
    echo        请先在 FEM_2d 项目根目录编译求解器
    pause
    exit /b 1
)

REM 检查 exe
if not exist "src-tauri\target\debug\femlab-studio.exe" (
    echo [错误] 未找到 femlab-studio.exe
    echo        请先运行: build_gui.bat
    pause
    exit /b 1
)

REM 复制 femcli 到应用目录(后端 locate_femcli 优先从 exe 同目录找)
copy /y "..\build\bin\femcli.exe" "src-tauri\target\debug\femcli.exe" >nul 2>&1

echo [1/1] 启动 FemLab Studio...
start "" "src-tauri\target\debug\femlab-studio.exe"

echo.
echo 已启动! 窗口标题: FemLab Studio - 二维杆系有限元分析
pause
