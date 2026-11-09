@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" goto dependencies
echo Przygotowanie DOSY Lab do pierwszego uruchomienia...
where py >nul 2>&1
if errorlevel 1 goto python_fallback
py -3 -m venv .venv
goto environment_check
:python_fallback
python -m venv .venv
:environment_check
if errorlevel 1 goto error
:dependencies
".venv\Scripts\python.exe" -c "import sys,numpy,scipy,PySide6,pyqtgraph,openpyxl; assert sys.version_info >= (3,11)" >nul 2>&1
if not errorlevel 1 goto launch
echo Instalowanie bibliotek. Internet jest potrzebny tylko do instalacji.
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto error
:launch
start "" ".venv\Scripts\pythonw.exe" "%~dp0app.py" %*
exit /b 0
:error
echo.
echo Nie udalo sie uruchomic DOSY Lab. Wymagany jest Python 3.11 lub nowszy.
pause
exit /b 1
