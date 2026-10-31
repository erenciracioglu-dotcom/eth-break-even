@echo off
rem Strata IQ3_S server on port 8080 for the ETH loop (ETH_Strata scheduled task). Output: logs\strata.log
cd /d D:\Strata
"D:\Strata\.venv\Scripts\python.exe" "D:\Strata\serve\server.py" --engine strata --config "D:\Strata\strata-iq3_s.json" --port 8080 --fit-max-tokens >> "D:\ETH-loop\logs\strata.log" 2>&1
