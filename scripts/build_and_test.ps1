# Build and run tests (PowerShell)
mkdir build
cd build
cmake ..
cmake --build . --config Debug
ctest --output-on-failure
