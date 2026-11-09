@echo off
setlocal
if exist "%~dp0.venv\Scripts\python.exe" goto dependencies
echo Setting up DOSY Lab for its first launch...
where py >nul 2>&1
if errorlevel 1 goto python_fallback
py -3 -m venv "%~dp0.venv"
goto environment_check
:python_fallback
python -m venv "%~dp0.venv"
:environment_check
if errorlevel 1 goto error
:dependencies
"%~dp0.venv\Scripts\python.exe" -c "import sys,numpy,scipy,PySide6,pyqtgraph,openpyxl; assert sys.version_info >= (3,11)" >nul 2>&1
if not errorlevel 1 goto launch
echo Installing dependencies. Internet access is only required for installation.
"%~dp0.venv\Scripts\python.exe" -m pip install -r "%~dp0bin\requirements.txt"
if errorlevel 1 goto error
:launch
start "" "%~dp0.venv\Scripts\pythonw.exe" "%~dp0bin\app.py" %*
exit /b 0
:error
echo.
echo DOSY Lab could not start. Python 3.11 or later is required.
pause
exit /b 1
