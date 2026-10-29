"""Final fixed-floor check and explicit disclosure of nonconverged peak gain."""
import json
from pathlib import Path
import numpy as np
from .reporting import refinement_note, outcome_interpretation, target_peak_ratios
from .authentication import authenticate_case, case_fingerprint
from d3d_optimization.case import (
    load_case, separatrix, pretrace, ShapeEvaluator, save_json, sha256)

import argparse
parser = argparse.ArgumentParser(description="Verify saved fixed-equilibrium D3D candidates")
parser.add_argument("--output",required=True,type=Path)
args = parser.parse_args()
p = args.output.resolve()
m = json.loads((p/"manifest.json").read_text(encoding="utf-8"))
wall = authenticate_case(p,m)
fingerprint = case_fingerprint(m)
cache = p/"finer_check.json"
d = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
if d.get("case_fingerprint") != fingerprint:
    assert sha256(m["shot_path"]) == m["shot_sha256"]
    wall = np.asarray(json.loads((p/"DIIID_geom.json").read_text())["limiter"])
    wall = np.vstack([wall, wall[0]])
    f,l,a = load_case(m["shot_path"],m["time_ms"])
    r,_,_ = separatrix(f,l,a)
    edges,paths,inside,_,drift = pretrace(f,r,a[1],wall,1280,.0005)
    e = ShapeEvaluator(f,l,wall,paths,inside,edges,r,a[1],m["target_original_segment"],641)
    result = e.evaluate(m["coefficients_m"])
    assert result["geometry_valid"]
    ratios = target_peak_ratios(
        [max(h.q[np.isin(h.segment,result["target_segments"])],default=0.) for h in result["maps"]],
        [max(b.q[np.isin(b.segment,e.baseline["target_segments"])],default=0.) for b in e.baseline["maps"]])
    if ratios is None:
        raise SystemExit("baseline_target_unilluminated: finest target peak comparison is undefined")
    d = {"name":"cubic_finer", "minimum_target_angle_deg":result["angle"],
         "feasible":result["feasible"], "target_peak_ratios":ratios,
         "all_wall_peak_ratios":[float(max(h.q)/b) for h,b in zip(result["maps"],e.baseline_peaks)],
         "max_path_psi_drift":drift, "unmapped_power":result["missing"],
         "angle_application":result["angle_application"], "constraint_penalties":result["constraint_penalties"]}
d.update(method="cubic",tubes=1280,step_m=.0005,surface_nodes=641,
         coefficients_m=m["coefficients_m"],shot_sha256=m["shot_sha256"],case_fingerprint=fingerprint)
assert d["feasible"] and d["minimum_target_angle_deg"] >= 1.5
assert max(d["all_wall_peak_ratios"]) <= 1+1e-10
save_json(cache,d)
checks = json.loads((p/"numerical_checks.json").read_text(encoding="utf-8"))
checks = [c for c in checks if c["name"] != "cubic_finer"]+[d]
save_json(p/"numerical_checks.json",checks)
m["optimization_outcome"].update(
    finest_target_peak_reduction_percent=100*(1-max(d["target_peak_ratios"])),
    target_gain_converged=False,
    interpretation=outcome_interpretation(checks))
m["search_script_sha256"] = sha256(Path(__file__).with_name("expanded_search.py"))
m["validation_script_sha256"] = sha256(Path(__file__).with_name("validate_expanded.py"))
m["final_refinement_script_sha256"] = sha256(Path(__file__))
save_json(p/"manifest.json",m)
report = p/"REPORT.md"
text = report.read_text(encoding="utf-8")
marker = "\n## 최종 세분화 확인과 개선율 해석\n"
text = text.split(marker)[0]
text += marker+"\n고정 후보를 cubic 1,280 tubes, 경로 간격 0.5 mm, 타깃 641 nodes로 한 단계 더 확인했습니다. 입사각 하한은 1.5°로 그대로 유지했습니다.\n\n"
text += f"최소 타깃 입사각 {d['minimum_target_angle_deg']:.6f}°, 전체 벽 peak ratio {max(d['all_wall_peak_ratios']):.6f}, shelf peak ratio {max(d['target_peak_ratios']):.6f}, shelf peak 감소 {100*(1-max(d['target_peak_ratios'])):.4f}%입니다.\n\n"
text += refinement_note(checks)
text += f"확대 탐색은 총 {m['counts']['evaluated']:,}개 고유 후보, 기하학적으로 유효한 {m['counts']['geometry_valid']:,}개 후보를 포함합니다. 재현 마지막 단계는 python -m d3d_optimization.final_refinement --output 결과폴더 입니다.\n"
report.write_text(text,encoding="utf-8")
print(m["optimization_outcome"])
