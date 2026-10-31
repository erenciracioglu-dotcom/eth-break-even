@echo off
rem Asks the ETH loop to stop after its current iteration; it then closes the Strata server it started.
echo stop> "D:\ETH-loop\STOP"
echo The ETH loop will stop after the current iteration.
timeout /t 5 >nul
