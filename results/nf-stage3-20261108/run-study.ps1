[CmdletBinding()]
param (
    [Parameter(Mandatory)] [string] $Study,
    [string] $Python = "python",
    [string] $GpuRuntime
)

$ErrorActionPreference = "Stop"
$studyPath = (Resolve-Path -LiteralPath $Study).Path
$repoPath = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "../..")).Path
$env:PYTHONPATH = (Resolve-Path -LiteralPath (Join-Path $repoPath ".local_dependencies/d3d-fixed-equilibrium")).Path
$env:PYTHONDONTWRITEBYTECODE = "1"
$env:PYTHONIOENCODING = "utf-8"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
if ($GpuRuntime) {
    $env:DIVERTOR_GPU_RUNTIME = (Resolve-Path -LiteralPath $GpuRuntime).Path
}
foreach ($pair in @(
    ,@("TEMP", "temp")
    ,@("TMP", "temp")
    ,@("MPLCONFIGDIR", "matplotlib")
    ,@("CUPY_CACHE_DIR", "cupy")
    ,@("CUDA_CACHE_PATH", "cuda")
)) {
    $cachePath = Join-Path $studyPath ("runtime/" + $pair[1])
    [void][IO.Directory]::CreateDirectory($cachePath)
    [Environment]::SetEnvironmentVariable($pair[0], $cachePath, "Process")
}

function Invoke-Stage ([string[]] $CliArgs) {
    Write-Output ("{0:u} {1}" -f [DateTime]::UtcNow, ($CliArgs -join " "))
    & $Python -B -u -m @CliArgs --study $studyPath
    if ($LASTEXITCODE -ne 0) {
        throw ("{0} exited with code {1}" -f ($CliArgs -join " "), $LASTEXITCODE)
    }
}

Push-Location -LiteralPath $repoPath
try {
    Invoke-Stage @("verification.stage3_campaign", "run")
    Invoke-Stage @("verification.stage3_campaign", "check")
    Invoke-Stage @("verification.stage3_report")
    Invoke-Stage @("verification.stage3_replay")
} finally {
    Pop-Location
}
