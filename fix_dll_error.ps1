# Windows DLL Error Auto-fix Script
# Purpose: Automatically install Visual C++ runtime and reinstall PySide6
# Usage: Run in PowerShell (right-click, select "Run with PowerShell")

Write-Host "====================================================" -ForegroundColor Cyan
Write-Host "  FEM_2d Windows DLL Error Auto-fix Tool" -ForegroundColor Cyan
Write-Host "====================================================" -ForegroundColor Cyan
Write-Host ""

# Check for admin privileges
$isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole] "Administrator")
if (-not $isAdmin) {
    Write-Host "[!] Note: This script needs admin rights to install Visual C++ runtime" -ForegroundColor Yellow
    Write-Host "[*] Recommended: Right-click this script and select 'Run with PowerShell'" -ForegroundColor Yellow
    Write-Host ""
}

Write-Host "[Step 1] Checking Visual C++ runtime..." -ForegroundColor Green

# Check if Visual C++ runtime is installed
$vcRedistPath = "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall"
$vcInstalled = $false

if (Test-Path $vcRedistPath) {
    $vcItems = Get-ChildItem $vcRedistPath | ForEach-Object { Get-ItemProperty $_.PSPath } | Where-Object { $_.DisplayName -like "*Visual C++*Redistributable*" }
    if ($vcItems) {
        Write-Host "[OK] Found Visual C++ runtime:" -ForegroundColor Green
        foreach ($item in $vcItems) {
            Write-Host "     - $($item.DisplayName)" -ForegroundColor Green
        }
        $vcInstalled = $true
    }
}

if (-not $vcInstalled) {
    Write-Host "[!] Visual C++ runtime not found - this is likely the cause of DLL errors" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "[Step 2] Downloading and installing Visual C++ runtime..." -ForegroundColor Green
    
    if ($isAdmin) {
        try {
            # Detect Python bitness
            $pythonPath = python -c "import sys; print(sys.executable)"
            $pythonBits = python -c "import struct; print(struct.calcsize('P') * 8)"
            
            Write-Host "[*] Detected Python: $pythonPath" -ForegroundColor Cyan
            Write-Host "[*] Python bitness: $pythonBits-bit" -ForegroundColor Cyan
            
            if ($pythonBits -eq 64) {
                $vcRedistUrl = "https://aka.ms/vs/17/release/vc_redist.x64.exe"
                $vcRedistFile = "$env:TEMP\vc_redist.x64.exe"
                $bits = "64"
            } else {
                $vcRedistUrl = "https://aka.ms/vs/17/release/vc_redist.x86.exe"
                $vcRedistFile = "$env:TEMP\vc_redist.x86.exe"
                $bits = "32"
            }
            
            Write-Host "[*] Downloading VC++ runtime ($bits-bit)..." -ForegroundColor Cyan
            $ProgressPreference = 'SilentlyContinue'
            Invoke-WebRequest -Uri $vcRedistUrl -OutFile $vcRedistFile -ErrorAction Stop
            
            Write-Host "[*] Installing VC++ runtime..." -ForegroundColor Cyan
            Start-Process -FilePath $vcRedistFile -ArgumentList "/install /quiet /norestart" -Wait -ErrorAction Stop
            
            Write-Host "[OK] VC++ runtime installation complete!" -ForegroundColor Green
            Remove-Item $vcRedistFile -Force -ErrorAction SilentlyContinue
        }
        catch {
            Write-Host "[!] Auto-download failed: $_" -ForegroundColor Red
            Write-Host "[*] Please manually download from: https://support.microsoft.com/en-us/help/2977003" -ForegroundColor Yellow
        }
    } else {
        Write-Host "[!] Admin rights required to install Visual C++ runtime" -ForegroundColor Red
        Write-Host "[*] Please right-click this script and select 'Run with PowerShell'" -ForegroundColor Yellow
    }
} else {
    Write-Host "[OK] Visual C++ runtime is installed, skipping installation" -ForegroundColor Green
}

Write-Host ""
Write-Host "[Step 3] Clearing pip cache and reinstalling PySide6..." -ForegroundColor Green

try {
    Write-Host "[*] Clearing pip cache..." -ForegroundColor Cyan
    python -m pip cache purge
    
    Write-Host "[*] Uninstalling old PySide6..." -ForegroundColor Cyan
    python -m pip uninstall -y PySide6
    
    Write-Host "[*] Installing latest PySide6..." -ForegroundColor Cyan
    python -m pip install --force-reinstall PySide6
    
    Write-Host "[OK] PySide6 reinstallation complete!" -ForegroundColor Green
}
catch {
    Write-Host "[!] Installation failed: $_" -ForegroundColor Red
}

Write-Host ""
Write-Host "[Step 4] Verifying fix..." -ForegroundColor Green

try {
    $output = python -c "from PySide6.QtWidgets import QApplication; print('OK')" 2>&1
    if ($output -eq "OK") {
        Write-Host "[OK] PySide6 import successful!" -ForegroundColor Green
        Write-Host ""
        Write-Host "====================================================" -ForegroundColor Green
        Write-Host "  Fix successful! You can now launch the GUI" -ForegroundColor Green
        Write-Host "====================================================" -ForegroundColor Green
        Write-Host ""
        Write-Host "Launch GUI: python quick_start_gui.py" -ForegroundColor Cyan
        Write-Host "Or run diagnostics: python diagnose.py" -ForegroundColor Cyan
    } else {
        Write-Host "[!] PySide6 import still has issues:" -ForegroundColor Yellow
        Write-Host $output -ForegroundColor Yellow
    }
}
catch {
    Write-Host "[!] Verification failed: $_" -ForegroundColor Red
}

Write-Host ""
Write-Host "[Press any key to exit...]" -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
