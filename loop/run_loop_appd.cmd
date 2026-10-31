@echo off
rem ETH loop driver for the Appendix D route, started by the ETH_Loop_AppD scheduled task. Output: logs\loop_appd_console.log
cd /d D:\ETH-loop
"D:\ETH-loop\.venv\Scripts\python.exe" -u "D:\ETH-loop\loop\loop_appd.py" >> "D:\ETH-loop\logs\loop_appd_console.log" 2>&1
