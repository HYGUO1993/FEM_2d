# build_exe.ps1
Write-Host "Building FEM_2d GUI Application..." -ForegroundColor Green

$ProjectRoot = $PSScriptRoot

if (Test-Path "$ProjectRoot\dist") { Remove-Item -Recurse -Force "$ProjectRoot\dist" }
if (Test-Path "$ProjectRoot\build\pyi") { Remove-Item -Recurse -Force "$ProjectRoot\build\pyi" }

$SolverPath = "$ProjectRoot\build\bin\Release\fem_run.exe"
if (-Not (Test-Path $SolverPath)) {
    Write-Host "fem_run.exe not found! Please build the C++ project in Release mode first." -ForegroundColor Red
    exit 1
}

Write-Host "Running PyInstaller..." -ForegroundColor Green
& "C:\Users\guoho\anaconda3\envs\fem_build_env\python.exe" -m PyInstaller --noconfirm `
    --collect-all PyQt6 `
    --name "FEM_2d_Studio" `
    --windowed `
    --onefile `
    --add-data "$SolverPath;bin" `
    --add-data "$ProjectRoot\example_continuous_beam.txt;." `
    --add-data "$ProjectRoot\example_multi_story_frame.txt;." `
    --add-data "$ProjectRoot\example_truss_bridge.txt;." `
    --workpath "$ProjectRoot\build\pyi" `
    --distpath "$ProjectRoot\dist" `
    "$ProjectRoot\python\fem_gui.py"

Write-Host "Packaging complete! Executable is located in dist/FEM_2d_Studio" -ForegroundColor Green
