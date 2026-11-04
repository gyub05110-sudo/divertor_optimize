param([switch]$NoOpen)
$ErrorActionPreference = 'Stop'
$viewerRoot = $PSScriptRoot
$viewerPort = 8780
$viewerUrl = "http://127.0.0.1:$viewerPort/"
$viewerPython = 'C:/Users/FPL-user/miniforge3/python.exe'
$alreadyServing = $false
try {
    $response = Invoke-WebRequest -Uri $viewerUrl -UseBasicParsing -TimeoutSec 2
    if ($response.Content -notmatch 'Divertor Shape Explorer') { throw 'Port 8780 belongs to another site.' }
    $alreadyServing = $true
} catch {
    if ($_.Exception.Message -eq 'Port 8780 belongs to another site.') { throw }
}
if (-not $alreadyServing) {
    $viewerProcess = Start-Process -FilePath $viewerPython -ArgumentList @(('"' + (Join-Path $viewerRoot 'serve.py') + '"')) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $viewerRoot 'server.stdout.log') -RedirectStandardError (Join-Path $viewerRoot 'server.stderr.log')
    @{pid=$viewerProcess.Id;started=$viewerProcess.StartTime.ToUniversalTime().ToString('o');url=$viewerUrl} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $viewerRoot 'server.json') -Encoding UTF8
}
Write-Output $viewerUrl
if (-not $NoOpen) { Start-Process $viewerUrl }
