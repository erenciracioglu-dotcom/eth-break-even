@echo off
rem Asks the Appendix D loop to stop after its current iteration; it then closes the Strata server it started.
echo stop> "D:\ETH-loop\STOP_APPD"
echo The Appendix D loop will stop after the current iteration.
timeout /t 5 >nul
