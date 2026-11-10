[CmdletBinding()]
param(
    [Parameter(Mandatory)] [string] $Study,
    [string] $Python = "python",
    [ValidateSet(1,2)] [int] $Workers = 2
)

$ErrorActionPreference = "Stop"
$root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "../..")).Path
$studyPath = (Resolve-Path -LiteralPath $Study).Path
$env:PYTHONPATH = (Resolve-Path -LiteralPath (Join-Path $root ".local_dependencies/d3d-fixed-equilibrium")).Path
$env:PYTHONDONTWRITEBYTECODE = "1"
$env:PYTHONIOENCODING = "utf-8"
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:TEMP = Join-Path $studyPath "runtime/temp"
[void][IO.Directory]::CreateDirectory($env:TEMP)
$env:TMP = $env:TEMP

function Invoke-Stage ([string[]] $CliArgs) {
    Write-Output ("{0:u} {1}" -f [DateTime]::UtcNow, ($CliArgs -join " "))
    & $Python -B -u -m @CliArgs --study $studyPath
    if ($LASTEXITCODE -ne 0) { throw ("Native campaign step failed ({0})" -f $LASTEXITCODE) }
}

Push-Location -LiteralPath $root
try {
    Invoke-Stage @("verification.stage5_dispatch", "run", "--workers", "$Workers")
    Invoke-Stage @("verification.stage5_dispatch", "check")
    Invoke-Stage @("verification.stage5_metrics")
} finally {
    Pop-Location
}
