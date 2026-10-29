$ErrorActionPreference='Stop'
$reviewPython=(Get-Command python).Source
$reviewArgs=@(('"'+(Join-Path $PSScriptRoot 'serve_review.py')+'"'),'--port','8777')
if (-not (Get-NetTCPConnection -LocalPort 8777 -State Listen -ErrorAction SilentlyContinue)) {
    Start-Process -FilePath $reviewPython -ArgumentList $reviewArgs -WindowStyle Hidden -RedirectStandardOutput (Join-Path $PSScriptRoot 'review-server.log') -RedirectStandardError (Join-Path $PSScriptRoot 'review-server-errors.log')
}
Write-Output 'http://127.0.0.1:8777/baseline_review/'
