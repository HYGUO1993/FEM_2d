@echo off
setlocal

set "ROOT=%~dp0.."
set "EXE=%ROOT%\build\bin\Release\fem_run.exe"
set "INPUT=%ROOT%\test05.txt"
set "OUTPUT=%ROOT%\build\Results.dat"

if not exist "%EXE%" (
  echo [ERROR] Release executable not found:
  echo %EXE%
  echo.
  echo Please build Release first:
  echo   cmake --build "%ROOT%\build" --config Release
  echo.
  pause
  exit /b 1
)

echo Running FEM_2d Release...
echo Input : %INPUT%
echo Output: %OUTPUT%
echo.
"%EXE%" --input "%INPUT%" --output "%OUTPUT%"
set "EXITCODE=%ERRORLEVEL%"
echo.
if "%EXITCODE%"=="0" (
  echo Done. Results written to:
  echo %OUTPUT%
) else (
  echo [ERROR] fem_run exited with code %EXITCODE%
)
echo.
pause
exit /b %EXITCODE%
