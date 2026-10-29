"""DIII-D fixed-equilibrium workflow integrated into the local FreeGS project."""
import argparse
from datetime import datetime
from pathlib import Path
import subprocess
import shutil
import sys
from uuid import uuid4

ROOT = Path(__file__).resolve().parent


def prepare_geometry(output):
    """Seed new outputs with the same pinned public geometry as the target tag."""
    output.mkdir(parents=True,exist_ok=True)
    for name in ("DIIID_geom.json","geometry_provenance.json"):
        target = output/name
        if target.exists():
            raise FileExistsError(f"geometry cache already exists: {target}")
        shutil.copyfile(ROOT/"d3d_optimization/data"/name,target)


def _run(module, output):
    subprocess.run([sys.executable, "-m", module, "--output", str(output)], cwd=ROOT, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="DIII-D fixed-equilibrium target-only heat-flux optimization")
    parser.add_argument("--mode", choices=["preflight", "optimize", "verify"], default="optimize")
    parser.add_argument("--design",choices=["floor-shelf","shelf"],default="floor-shelf",
                        help="default: latest eight-variable floor+shelf design; shelf: historical four-variable design")
    parser.add_argument("--shot", type=Path, default=ROOT/"input/D3D/d3d/d3d_shot_1a3155370c.parquet")
    parser.add_argument("--target-config", type=Path, default=ROOT/"d3d_optimization/target_reference.json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--time-ms", type=float, default=1400.)
    parser.add_argument("--bound-mm", type=float, default=100.)
    parser.add_argument("--min-cases", type=int, default=1000)
    parser.add_argument("--generations", type=int, default=60)
    args = parser.parse_args(argv)
    if args.mode == "verify" and args.output is None:
        parser.error("--mode verify requires --output pointing to an existing result")
    if args.min_cases < 1000 or args.generations < 0 or not 0 < args.bound_mm <= 1000:
        parser.error("min-cases >=1000, generations >=0 and 0<bound-mm<=1000 required")
    output = args.output.resolve() if args.output else ROOT/"input/D3D/optimized"/(
        "run-"+datetime.now().strftime("%Y%m%d-%H%M%S")+"-"+uuid4().hex[:8])
    if args.mode != "verify" and output.exists():
        parser.error("new runs require a new output path; existing results are preserved")
    from d3d_optimization import case
    if args.mode == "verify":
        import json
        manifest = json.loads((output/"manifest.json").read_text(encoding="utf-8"))
        if "config" in manifest and "surfaces" in manifest["config"]:
            if case.sha256(args.shot) != manifest["shot_sha256"]:
                parser.error("verification shot hash differs from saved case")
            if manifest.get("preflight_only"):
                parser.error("this is a baseline preflight, not a searched candidate")
            _run("d3d_optimization.multisurface_validate",output)
            return 0
        if "counts" in manifest:
            if case.sha256(args.shot) != manifest["shot_sha256"]:
                parser.error("verification shot hash differs from saved case")
            config = json.loads(args.target_config.read_text(encoding="utf-8"))
            if config != manifest["target_configuration"]:
                parser.error("verification target config differs from saved case")
            _run("d3d_optimization.validate_expanded", output)
            _run("d3d_optimization.final_refinement", output)
        else:
            case.main(["--shot",str(args.shot),"--output",str(output),"--target-config",str(args.target_config),"--verify-only"])
        return 0
    common = ["--shot",str(args.shot),"--target-config",str(args.target_config),"--time-ms",str(args.time_ms),"--preflight-only"]
    if args.design == "floor-shelf":
        if args.target_config.resolve() != (ROOT/"d3d_optimization/target_reference.json").resolve():
            parser.error("custom --target-config requires --design shelf")
        output.mkdir(parents=True)
        reference=output/"reference_case"
        prepare_geometry(reference)
        case.main(common+["--output",str(reference)])
        command=[sys.executable,"-m","d3d_optimization.multisurface_search","--output",str(output),
                 "--reference-case",str(reference),"--min-cases",str(args.min_cases),
                 "--generations",str(args.generations),"--bound-mm",str(args.bound_mm)]
        if args.mode == "preflight": command.append("--preflight-only")
        subprocess.run(command,cwd=ROOT,check=True)
        if args.mode == "optimize": _run("d3d_optimization.multisurface_validate",output)
        print(f"Floor+shelf research result: {output}")
        return 0
    if args.mode == "preflight":
        prepare_geometry(output)
        case.main(common+["--output",str(output)])
        return 0
    output.mkdir(parents=True)
    reference = output/"reference_case"
    prepare_geometry(reference)
    case.main(common+["--output",str(reference)])
    import json
    reference_manifest = json.loads((reference/"manifest.json").read_text(encoding="utf-8"))
    if reference_manifest["angle_application"]["status"] == "target_unilluminated":
        print(f"target_unilluminated: no powered first hits on the configured baseline target; {reference}")
        return 2
    subprocess.run([sys.executable,"-m","d3d_optimization.expanded_search",
                    "--output",str(output),"--reference-case",str(reference),
                    "--min-cases",str(args.min_cases),"--generations",str(args.generations),
                    "--bound-mm",str(args.bound_mm)],cwd=ROOT,check=True)
    _run("d3d_optimization.validate_expanded",output)
    _run("d3d_optimization.final_refinement",output)
    print(f"Verified research result: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
