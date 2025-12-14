# Build and run tests (PowerShell)
mkdir build
cd build
$cmake = "${env:ProgramFiles}\CMake\bin\cmake.exe"
if (-not (Test-Path $cmake)) { $cmake = "${env:ProgramFiles(x86)}\CMake\bin\cmake.exe" }
$ctest = Join-Path (Split-Path $cmake) 'ctest.exe'
& $cmake ..
& $cmake --build . --config Debug
& $ctest -C Debug --output-on-failure
