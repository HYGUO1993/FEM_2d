@echo off
setlocal

set "ROOT=%~dp0.."
set "RESULTS=%ROOT%\build\Results.dat"
set "OUT=%ROOT%\build\plot.png"
set "PY_CMD="

where python >nul 2>nul
if %ERRORLEVEL%==0 (
  set "PY_CMD=python"
) else (
  where py >nul 2>nul
  if %ERRORLEVEL%==0 (
    set "PY_CMD=py -3"
  )
)

if "%PY_CMD%"=="" (
  echo [ERROR] Python not found. Please install Python 3 and try again.
  pause
  exit /b 1
)

if not exist "%RESULTS%" (
  echo [INFO] Results file not found. Running solver first...
  call "%~dp0run_release.bat"
  if errorlevel 1 (
    echo [ERROR] Failed to generate results file.
    pause
    exit /b 1
  )
)

echo Generating visualization...
call %PY_CMD% "%~dp0visualize_results.py" --results "%RESULTS%" --scale 1000 --out "%OUT%"
if errorlevel 1 (
  echo [ERROR] Visualization failed. If modules are missing, install:
  echo   pip install matplotlib numpy
  pause
  exit /b 1
)

echo Done. Image saved to:
echo %OUT%
start "" "%OUT%"
echo.
pause
exit /b 0
