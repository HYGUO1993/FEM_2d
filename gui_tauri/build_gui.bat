@echo off
REM ============================================================
REM FemLab Studio (Tauri) - 构建脚本
REM 编译 Rust 后端, 并复制 femcli 到产物目录
REM 打包发布: 见 README 或 tauri build
REM ============================================================
setlocal
cd /d "%~dp0"

echo ============================================================
echo  FemLab Studio - 构建 (cargo build)
echo ============================================================

REM 检查 femcli
if not exist "..\build\bin\femcli.exe" (
    echo [警告] 未找到 ..\build\bin\femcli.exe
    echo        请先在 FEM_2d 项目根目录编译求解器
)

echo [1/2] cargo build...
cd src-tauri
call cargo build
if errorlevel 1 (
    echo [错误] cargo build 失败
    pause
    exit /b 1
)
cd ..

REM 复制 femcli 到 debug 产物目录
if exist "..\build\bin\femcli.exe" (
    copy /y "..\build\bin\femcli.exe" "src-tauri\target\debug\femcli.exe" >nul
    echo [2/2] femcli.exe 已复制到 target\debug\
)

echo.
echo 构建完成! 运行 dev.bat 启动开发版
pause
