[CmdletBinding()]
param(
    [string] $Stage3 = "D:/DivertorOptimization/studies/nf-stage3-20261108",
    [string] $OutputDirectory = $PSScriptRoot,
    [string] $Python = "python"
)

$ErrorActionPreference = "Stop"
$repoPath = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "../../..")).Path
$oldPythonPath = $env:PYTHONPATH
Push-Location -LiteralPath $repoPath
try {
    $env:PYTHONPATH = (Resolve-Path -LiteralPath ".local_dependencies/d3d-fixed-equilibrium").Path
    $xOutput = Join-Path $OutputDirectory "lower-x-diagnostics.json"
    & $Python -B -m verification.stage4_xpoint_diagnostics --stage3 $Stage3 --output $xOutput
    if ($LASTEXITCODE -ne 0) { throw "Lower-X development diagnosis failed" }

    $wallOutput = Join-Path $OutputDirectory "fixed-wall-regression.json"
    @'
from hashlib import sha256
import json
from pathlib import Path
import sys

from d3d_optimization.user_input import load_user_input
from d3d_optimization.user_precheck import fixed_wall_precheck

case_id = "b5e07dde39"
source = Path(sys.argv[1])/"inputs"/"raw"/f"d3d_shot_{case_id}.parquet"
expected = "a783f1ae24d01a159ad1123e24dd51f99ba6f60dd3f7577851b28d1229b0b749"
if sha256(source.read_bytes()).hexdigest() != expected:
    raise ValueError("archived development source changed")
user = load_user_input(source, 1400.)
check = fixed_wall_precheck(user)
if not check["excluded"] or [r["original_material_segment"] for r in check["fixed_segment_intersections"]] != [4, 5]:
    raise ValueError("fixed material exclusion failed to reproduce")
report = dict(case_id=case_id, equilibrium_sha256=expected, wall_sha256=user.provenance["wall_sha256"],
              time_ms=user.time_ms, fixed_wall=check)
Path(sys.argv[2]).write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf8")
print(case_id, check["classification"], "original material segments 4 and 5")
'@ | & $Python -B - $Stage3 $wallOutput
    if ($LASTEXITCODE -ne 0) { throw "Fixed-wall regression failed" }
} finally {
    $env:PYTHONPATH = $oldPythonPath
    Pop-Location
}
