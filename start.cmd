@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\start.ps1" %*
if errorlevel 1 (
  echo.
  echo Startup failed. See the message above and logs\launcher\ for details.
  pause
  exit /b 1
)
