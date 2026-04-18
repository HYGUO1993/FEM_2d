@echo off
REM Windows DLL Error Auto Fix Script Wrapper
REM Purpose: Start PowerShell script with proper permissions
REM Usage: Double-click this .bat file to run

setlocal enabledelayedexpansion

REM Get current script directory
set SCRIPT_DIR=%~dp0

REM Check if fix_dll_error.ps1 exists
if not exist "%SCRIPT_DIR%fix_dll_error.ps1" (
    echo [ERROR] Cannot find fix_dll_error.ps1
    echo Please ensure this file is in the same directory as fix_dll_error.ps1
    pause
    exit /b 1
)

REM Run PowerShell script (with elevated privileges)
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '%SCRIPT_DIR%fix_dll_error.ps1'"

pause
