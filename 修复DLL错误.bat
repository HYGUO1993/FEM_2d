@echo off
REM Windows DLL 错误自动修复脚本包装
REM 用途: 以正确的权限启动 PowerShell 脚本
REM 使用方法: 直接双击此 .bat 文件运行

setlocal enabledelayedexpansion

REM 获取当前脚本目录
set SCRIPT_DIR=%~dp0

REM 检查 fix_dll_error.ps1 是否存在
if not exist "%SCRIPT_DIR%fix_dll_error.ps1" (
    echo [ERROR] 找不到 fix_dll_error.ps1
    echo 请确保此文件与 fix_dll_error.ps1 在同一目录
    pause
    exit /b 1
)

REM 运行 PowerShell 脚本（提升权限）
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%SCRIPT_DIR%fix_dll_error.ps1'"

pause
