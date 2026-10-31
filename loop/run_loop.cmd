@echo off
rem ETH loop driver, started by the ETH_Loop scheduled task. Console output: logs\loop_console.log
cd /d D:\ETH-loop
"D:\ETH-loop\.venv\Scripts\python.exe" -u "D:\ETH-loop\loop\loop.py" >> "D:\ETH-loop\logs\loop_console.log" 2>&1
