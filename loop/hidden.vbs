' Runs the .cmd file given as the first argument without a visible window and waits for it.
' Used by the ETH_Strata and ETH_Loop scheduled tasks.
Set sh = CreateObject("WScript.Shell")
rc = sh.Run("cmd.exe /c """ & WScript.Arguments(0) & """", 0, True)
WScript.Quit rc
