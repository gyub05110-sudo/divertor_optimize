"""Select a candidate passing the same 1.5 deg floor in all numeric checks."""
import csv
import json
from pathlib import Path
import numpy as np
from .authentication import authenticate_case, case_fingerprint
from .reporting import target_peak_ratios
from d3d_optimization.case import (
    load_case, separatrix, pretrace, ShapeEvaluator, rows_for, plots,
    write_angle_scope, save_json, WIDTHS_MM, sha256)

import argparse
parser = argparse.ArgumentParser(description="Verify saved fixed-equilibrium D3D candidates")
parser.add_argument("--output",required=True,type=Path)
args = parser.parse_args()
p = args.output.resolve()
m = json.loads((p/"manifest.json").read_text(encoding="utf-8"))
wall = authenticate_case(p,m)
evaluators = []
for name, method, tubes, step, nodes in [
        ("cubic_coarse", "cubic", 320, .002, 161),
        ("cubic_refined", "cubic", 640, .001, 321),
        ("quintic_refined", "quintic", 640, .001, 321),
        ("cubic_finer", "cubic", 1280, .0005, 641)]:
    field, lcfs, axis = load_case(m["shot_path"], m["time_ms"], method)
    r_sep, xp, _ = separatrix(field, lcfs, axis)
    edges, paths, inside, _, drift = pretrace(field, r_sep, axis[1], wall, tubes, step)
    e = ShapeEvaluator(field, lcfs, wall, paths, inside, edges, r_sep, axis[1], m["target_original_segment"], nodes)
    evaluators.append((name, e, drift))

records = [json.loads(line) for line in (p/"evaluations.jsonl").read_text(encoding="utf-8").splitlines()]
assert len({tuple(r["coefficients_m"]) for r in records}) == len(records)
assert sum(r["geometry_valid"] for r in records) >= 1000
ranked = sorted((r for r in records if r["feasible"]), key=lambda r: (
    r["worst_all_wall_peak_ratio"], r["worst_target_peak_ratio"], r["mean_target_peak_ratio"]))
attempts = []
chosen = None
for rank, row in enumerate(ranked):
    checks = []
    results = []
    for name, e, drift in evaluators:
        result = e.evaluate(row["coefficients_m"])
        if not result["geometry_valid"]:
            checks.append({"name": name, "geometry_valid": False, "feasible": False})
            break
        peaks = [float(max(h.q)) for h in result["maps"]]
        ratios = [float(q/b) for q, b in zip(peaks, e.baseline_peaks)]
        target_peaks = [float(max(h.q[np.isin(h.segment, result["target_segments"])], default=0.)) for h in result["maps"]]
        base_target = [float(max(h.q[np.isin(h.segment, e.baseline["target_segments"])], default=0.)) for h in e.baseline["maps"]]
        target_ratios = target_peak_ratios(target_peaks,base_target)
        if target_ratios is None:
            checks.append({"name":name,"geometry_valid":True,"feasible":False,
                           "status":"baseline_target_unilluminated","target_peak_ratios":None,
                           "all_wall_peak_ratios":ratios})
            break
        checks.append({"name": name, "method": e.field.method, "tubes": len(e.paths),
                       "minimum_target_angle_deg": result["angle"], "feasible": result["feasible"],
                       "all_wall_peak_ratios": ratios, "target_peak_ratios": target_ratios,
                       "unmapped_power": result["missing"], "max_path_psi_drift": drift,
                       "angle_application": result["angle_application"],
                       "constraint_penalties": result["constraint_penalties"]})
        results.append(result)
        # Non-worsening all-wall heat is required in every checked representation.
        if not result["feasible"] or max(ratios) > 1+1e-10:
            break
    passed = len(checks) == len(evaluators) and all(c["feasible"] and max(c["all_wall_peak_ratios"]) <= 1+1e-10 for c in checks)
    attempts.append({"case": row["case"], "rank": rank, "passed": passed, "checks": checks})
    if rank % 64 == 0:
        print(f"verification rank {rank}; passed={passed}", flush=True)
    if passed:
        chosen = row, checks, results[0]
        break
if chosen is None:
    m["numerically_verified_feasible"] = False
    save_json(p/"manifest.json",m)
    save_json(p/"candidate_verification.json",attempts)
    raise SystemExit("No candidate passed the fixed 1.5 degree constraint in every numerical check.")
row, checks, result = chosen
e = evaluators[0][1]
rows = rows_for(e, e.baseline, "baseline")+rows_for(e, result, "candidate")
with (p/"width_scan.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
np.savez_compressed(p/"case_arrays.npz", R=e.field.r, Z=e.field.z, psi=e.field.psi, lcfs=e.lcfs,
                    baseline_wall=wall, candidate_wall=result["wall"], target=result["target"],
                    coefficients=row["coefficients_m"], x_sources=e.sources[0].x)
field, lcfs, axis = load_case(m["shot_path"], m["time_ms"])
_, xp, _ = separatrix(field, lcfs, axis)
plots(p, field, lcfs, xp, e, result)
write_angle_scope(p, e, result)
m.update(coefficients_m=row["coefficients_m"], results=rows,
         minimum_target_angle_deg=result["angle"], angle_application=result["angle_application"],
         constraint_penalties=result["constraint_penalties"], numerically_verified_feasible=True,
         coarse_best_replaced_after_verification=(row["coefficients_m"] != m["coefficients_m"]),
         retained_case=row["case"], verification_attempt_count=len(attempts),
         max_actual_displacement_m=row["max_displacement_m"],
         optimization_outcome={"worst_peak_ratio": row["worst_all_wall_peak_ratio"],
                               "worst_target_peak_ratio": row["worst_target_peak_ratio"],
                               "mean_target_peak_ratio": row["mean_target_peak_ratio"],
                               "peak_reduction_percent": 100*(1-row["worst_all_wall_peak_ratio"]),
                               "worst_target_peak_reduction_percent": 100*(1-row["worst_target_peak_ratio"])})
save_json(p/"manifest.json", m)
save_json(p/"numerical_checks.json", checks)
finer = dict(checks[-1],coefficients_m=m["coefficients_m"],shot_sha256=m["shot_sha256"],case_fingerprint=case_fingerprint(m))
save_json(p/"finer_check.json",finer)
save_json(p/"candidate_verification.json", attempts)
counts = m["counts"]
lines = ["# DIII-D 확대 범위 형상 최적화", "",
         f"서로 다른 {counts['evaluated']:,}개 후보를 평가했습니다. 기하학적으로 유효한 후보 {counts['geometry_valid']:,}개, 초기 계산에서 고정 제약을 만족한 후보 {counts['feasible']:,}개입니다. 요구한 최소 1,000개는 기하학적으로 유효한 후보 수로 충족했습니다.", "",
         "## 범위와 목적", "",
         f"평형은 {m['time_ms']:g} ms 원본 EFIT로 고정했습니다. 명시적 타깃 {m['target_configuration']['name']} (원본 segment {m['target_original_segment']})의 끝점·끝점 기울기는 고정하고 계수 상한 ±{1000*m['coefficient_bound_m']:g} mm를 사용했습니다. 초기 공간 샘플링은 세 단계에서 각각 1,024개이며 세 초기값으로 differential evolution을 추가 실행했습니다. 길이 증가는 기준 대비 20% 이내이며 벽·LCFS·자기 교차를 금지합니다.", "",
         "1.5° 하한은 지정 타깃의 양의 파워 첫 충돌점에만 적용합니다. 비타깃 벽, 차폐된 후속 교차점, 파워 0인 점은 최소각과 각도 벌점에서 제외합니다. 폭은 1·2·3·5·10 mm, 유입 outer branch 파워는 1입니다.", "",
         "기본 목적은 폭 전체의 최대 전체 벽 peak ratio 최소화입니다. 이 값이 같은 후보들에서는 최대 shelf peak ratio, 이어 평균 shelf peak ratio가 낮은 후보를 선택했습니다. shelf 개선을 전체 벽 개선으로 해석하지 않습니다. 탐색용 scalar 점수와 최종 후보의 사전식 선택 규칙은 구분하며, 최종 후보는 모든 수치 확인에서 전체 벽 최대값이 증가하지 않아야 합니다.", "",
         "## 유지한 후보", "",
         f"원본 평가 case {row['case']}; spline 계수 (m): {row['coefficients_m']}. 실제 최대 표면 변위는 {1000*row['max_displacement_m']:.3f} mm입니다. manifest의 coefficient_bound_m는 계수 탐색 상한이며 실제 선택 변위를 뜻하지 않습니다.", "",
         f"전체 벽 최대 열유속 감소: {100*(1-row['worst_all_wall_peak_ratio']):.4f}%. shelf의 worst-case peak ratio 감소: {100*(1-row['worst_target_peak_ratio']):.4f}% (초기 계산 기준).", "",
         "| λq (mm) | 기준 shelf q/P (m⁻²) | 후보 shelf q/P (m⁻²) | shelf 감소 (%) | 전체 벽 후보/기준 | 후보 타깃 파워 |", "|---:|---:|---:|---:|---:|---:|"]
for b,c in zip(rows[:5], rows[5:]):
    lines.append(f"| {c['lambda_q_mm']} | {b['target_peak_m2inv']:.7g} | {c['target_peak_m2inv']:.7g} | {100*(1-c['target_peak_m2inv']/b['target_peak_m2inv']):.4f} | {c['all_wall_peak_ratio']:.6f} | {c['target_power']:.6g} |")
lines += ["", "## 같은 1.5° 하한에서 수치 확인", "",
         f"초기 순위에서 {len(attempts)}개 후보를 확인해 cubic 320/640/1280 tubes 및 quintic 640 tubes 모두에서 제약을 통과하는 첫 후보를 유지했습니다. 하한값을 바꾸는 민감도 스캔은 수행하지 않았습니다.", "",
          "| 표현 | 최소 타깃 입사각 (°) | 전체 벽 worst ratio | shelf worst ratio | 미포착 파워 |", "|---|---:|---:|---:|---:|"]
for c in checks:
    lines.append(f"| {c['name']} | {c['minimum_target_angle_deg']:.6f} | {max(c['all_wall_peak_ratios']):.6f} | {max(c['target_peak_ratios']):.6f} | {c['unmapped_power']:.3g} |")
lines += ["", "## 결과의 한계", "",
          "원본 ψ 격자는 약 26.6×50 mm로, mm 단위 형상·열유속 수치가 실제 장치의 해당 정확도를 의미하지 않습니다. 공개 형상과 실제 샷 캠페인의 대응 및 하드웨어 정체는 미확인입니다. 비타깃 벽 열부하를 계속 포함했고 포착된 파워를 재정규화하지 않았습니다. 실제 MW/m², detachment·ELM·재료 온도·제작 가능성을 예측하지 않습니다.", "",
          "전체 벽 최대값이 개선되지 않았다면 이는 전체 열부하 개선 형상이 아닙니다. 열유속의 독립적인 표면 적분 검증 및 젖음 면적 지표가 아직 없어 인증된 열 최적값으로 주장하지 않습니다.", "",
          "재현: `python optimize_d3d.py --mode optimize` 및 `--mode verify`. 모든 후보: `evaluations.jsonl`; 입력 hash·설정: `manifest.json`; 후보 확인: `candidate_verification.json`; 점별 제약 범위: `angle_scope.csv`."]
(p/"REPORT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
print("VERIFIED", m["optimization_outcome"], "angles", [c["minimum_target_angle_deg"] for c in checks], flush=True)
