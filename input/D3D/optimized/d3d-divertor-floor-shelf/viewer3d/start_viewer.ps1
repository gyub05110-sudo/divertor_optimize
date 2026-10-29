$ErrorActionPreference = 'Stop'
$viewerDirectory = $PSScriptRoot
$viewerPort = 8776
$viewerUrl = "http://127.0.0.1:$viewerPort/"
$viewerListener = Get-NetTCPConnection -LocalPort $viewerPort -State Listen -ErrorAction SilentlyContinue
if (-not $viewerListener) {
    $viewerPython = (Get-Command python).Source
    $viewerArguments = @(('"' + (Join-Path $viewerDirectory 'server.py') + '"'), '--port', "$viewerPort")
    Start-Process -FilePath $viewerPython -ArgumentList $viewerArguments -WindowStyle Hidden -RedirectStandardOutput (Join-Path $viewerDirectory 'server.log') -RedirectStandardError (Join-Path $viewerDirectory 'server-errors.log')
    Start-Sleep -Milliseconds 800
}
$viewerExpected = Get-Content -LiteralPath (Join-Path $viewerDirectory 'data.json') -Raw | ConvertFrom-Json
$viewerServed = Invoke-RestMethod -Uri ($viewerUrl + 'data.json') -TimeoutSec 5
if ($viewerServed.source_geometry_sha256 -ne $viewerExpected.source_geometry_sha256 -or $viewerServed.retained_case -ne $viewerExpected.retained_case) {
    throw 'Port 8776 is serving a different page. Close that service or use another port before launching this viewer.'
}
$viewerChrome = 'C:\Program Files\Google\Chrome\Application\chrome.exe'
Start-Process -FilePath $viewerChrome -ArgumentList $viewerUrl
Write-Output $viewerUrl
