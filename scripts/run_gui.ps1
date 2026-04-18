#!/usr/bin/env pwsh
<#
.SYNOPSIS
    FEM_2d GUI Application Launcher (PowerShell)

.DESCRIPTION
    This script launches the FEM_2d GUI application.
    It checks for required dependencies and handles environment setup.

.EXAMPLE
    .\scripts\run_gui.ps1

.NOTES
    Requirements:
      - Python 3.8+ with conda or system Python
      - PySide6, matplotlib, numpy packages
#>

param(
    [switch]$SkipCheck = $false
)

# Get the project root directory
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$PythonPath = Join-Path $ProjectRoot "python"

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "FEM_2d GUI Application Launcher" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""

# Check if Python is available
try {
    $PythonVersion = python --version 2>&1
    Write-Host "[INFO] Python version: $PythonVersion"
} catch {
    Write-Host "[ERROR] Python not found in PATH" -ForegroundColor Red
    Write-Host "Please install Python 3.8+ or activate a conda environment with Python" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host "[INFO] Project root: $ProjectRoot"
Write-Host ""

if (-not $SkipCheck) {
    Write-Host "[INFO] Checking dependencies..."
    
    $MissingPackages = @()
    
    # Check PySide6
    try {
        python -c "import PySide6" 2>$null
    } catch {
        $MissingPackages += "PySide6"
    }
    
    # Check matplotlib
    try {
        python -c "import matplotlib" 2>$null
    } catch {
        $MissingPackages += "matplotlib"
    }
    
    # Check numpy
    try {
        python -c "import numpy" 2>$null
    } catch {
        $MissingPackages += "numpy"
    }
    
    if ($MissingPackages.Count -gt 0) {
        Write-Host "[WARNING] Missing packages: $($MissingPackages -join ', ')" -ForegroundColor Yellow
        Write-Host ""
        Write-Host "[INFO] Installing missing packages..." -ForegroundColor Cyan
        
        & pip install -q --no-warn-script-location @MissingPackages
        
        if ($LASTEXITCODE -ne 0) {
            Write-Host "[ERROR] Failed to install dependencies" -ForegroundColor Red
            Write-Host "Please install manually with:" -ForegroundColor Yellow
            Write-Host "  pip install PySide6 matplotlib numpy" -ForegroundColor Gray
            Read-Host "Press Enter to exit"
            exit 1
        }
    }
    
    Write-Host "[INFO] All dependencies available" -ForegroundColor Green
    Write-Host ""
}

# Check if FEM release executable is built
$FemExe = Join-Path $ProjectRoot "build" "bin" "Release" "fem_run.exe"
if (-not (Test-Path $FemExe)) {
    Write-Host "[WARNING] Release executable not found" -ForegroundColor Yellow
    Write-Host "[INFO] The solver will not work until you build the Release version:" -ForegroundColor Cyan
    Write-Host "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release" -ForegroundColor Gray
    Write-Host "  cmake --build build --config Release" -ForegroundColor Gray
    Write-Host ""
}

Write-Host "[INFO] Launching GUI..." -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""

# Launch the GUI
Set-Location $ProjectRoot
& python "$PythonPath\fem_gui.py"

Write-Host ""
Write-Host "[INFO] GUI closed" -ForegroundColor Green
