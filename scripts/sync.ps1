# PowerShell wrapper: scripts\sync.ps1 push|pull|check  (delegates to Git Bash sync.sh)
param([Parameter(Mandatory=$true)][ValidateSet('push','pull','check')][string]$Action)
$bash = "C:\Program Files\Git\bin\bash.exe"
if (-not (Test-Path $bash)) { $bash = (Get-Command bash).Source }
& $bash "$PSScriptRoot/sync.sh" $Action
exit $LASTEXITCODE
