param(
    [Parameter(Mandatory)][string]$Python,
    [Parameter(Mandatory)][string]$OtherEquilibrium,
    [Parameter(Mandatory)][string]$Study,
    [ValidateSet("Search","Replay","Analyze")][string]$Phase = "Search"
)
$ErrorActionPreference = "Stop"
$repoPath = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$studyPath = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Study)
$otherPath = (Resolve-Path -LiteralPath $OtherEquilibrium).Path
$pairs = @(
    @{name="e806f47778"; path=$otherPath; sha="dc0eba5720dc0a8b9894ae6fdaf5f346b4d5b2a06ca4d40cd83be8def591f0a5"; order=@("full","floor-only")},
    @{name="1a3155370c"; path=(Join-Path $repoPath "input/D3D/d3d/d3d_shot_1a3155370c.parquet"); sha="cae74032203593dfc1e8d9b7571785a3b3ad558f6e4265e938a33f16e422022c"; order=@("floor-only","full")}
)
Push-Location -LiteralPath $repoPath
try {
    foreach ($pair in $pairs) {
        if ((Get-FileHash -LiteralPath $pair.path -Algorithm SHA256).Hash -ne $pair.sha) {
            throw "Incorrect equilibrium bytes for $($pair.name)"
        }
        $caseDir = Join-Path $studyPath $pair.name
        $full = Join-Path $caseDir "full"
        $floor = Join-Path $caseDir "floor-only"
        $replay = Join-Path $caseDir "replay"
        if ($Phase -eq "Search") {
            foreach ($family in $pair.order) {
                $out = Join-Path $caseDir $family
                & $Python -B optimize_divertor.py --equilibrium $pair.path --time-ms 1400 --backend gpu --budget-seconds 360 --seed 20711 --family $family --output $out
                if ($LASTEXITCODE -ne 0) { throw "Search failed in $out" }
            }
        } elseif ($Phase -eq "Replay") {
            & $Python -B -m verification.user_shape_replay --full $full --floor $floor --equilibrium $pair.path --output $replay --backend gpu --levels cubic1280 quintic1280 cubic640
            if ($LASTEXITCODE -ne 0) { throw "Replay failed in $replay" }
        } else {
            & $Python -B -m verification.user_shape_results --full $full --floor $floor --replay $replay --output (Join-Path $caseDir "analysis")
            if ($LASTEXITCODE -ne 0) { throw "Analysis failed in $caseDir" }
        }
    }
} finally {
    Pop-Location
}
