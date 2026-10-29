"""Reproducible research case: fixed DIII-D EFIT, one angle floor, width scan.

Run with python -m divertor_optimization.d3d_frozen_case --shot PATH --output DIR --target-config JSON.
PyArrow is an optional dependency required only for this parquet case runner.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import time
import urllib.request

import numpy as np
from scipy.optimize import brentq, root, differential_evolution
import shapely
from shapely.geometry import LineString, Polygon

from .frozen_field import FrozenField, exponential_source
from .frozen_mapping import (Hits, first_hits, map_heat, deform_target, trace_paths,
                             target_angle_status, ANGLE_SCOPE)

ANGLE_FLOOR = 1.5
WIDTHS_MM = np.array([1, 2, 3, 5, 10])
GEOMETRY_PATH = "src/examples/TokaMaker/DIIID/DIIID_geom.json"


def candidate_feasible(angle, missing, power_error):
    return bool(angle is not None and np.isfinite(angle) and angle >= ANGLE_FLOOR and
                np.isfinite(missing) and 0 <= missing <= 1e-3 and
                np.isfinite(power_error) and 0 <= power_error <= 1e-4)


def select_target(wall, geometry_hash, config):
    """Require an explicit, coordinate-checked target; never infer it from heat."""
    if not isinstance(config, dict) or not config.get("name"):
        raise ValueError("explicit named target configuration is required")
    if config.get("geometry_sha256") != geometry_hash:
        raise ValueError("target configuration geometry hash does not match")
    index = config.get("segment_index")
    if type(index) is not int or not 0 < index < len(wall)-2:
        raise ValueError("target segment must be an interior original contour segment")
    expected = np.asarray(config.get("endpoints_RZ_m"), dtype=float)
    if expected.shape != (2, 2) or not np.allclose(expected, wall[index:index+2], rtol=0, atol=1e-8):
        raise ValueError("target segment endpoints do not match the public contour")
    if type(config.get("physical_identity_validated")) is not bool:
        raise ValueError("target physical identity validation status must be explicit")
    return index


def constraint_penalties(angle_status, missing):
    """Unilluminated target is an unusable trial, never an angle violation."""
    angle = angle_status["minimum_angle_deg"]
    return {"angle_lower_bound": 100*max(0, ANGLE_FLOOR-angle) if angle is not None else 0.,
            "target_unilluminated": 1000. if angle_status["status"] == "target_unilluminated" else 0.,
            "invalid_target_angle": 1000. if angle_status["status"] == "invalid_target_angle" else 0.,
            "unmapped_power": 1000*max(0, missing-1e-3)}


def angle_text(angle, status=None):
    if angle is not None:
        return f"{angle:.4f}"
    return "N/A (invalid target angle)" if status == "invalid_target_angle" else "N/A (no powered target first hits)"


def case_status(result):
    if result["feasible"]:
        return "FEASIBLE"
    status = result["angle_application"]["status"]
    return status.upper() if status in ("target_unilluminated", "invalid_target_angle") else "INFEASIBLE UNDER FIXED CONSTRAINTS"


def verified_feasibility(coarse_feasible, checks):
    return bool(coarse_feasible and checks and all(c.get("feasible", False) for c in checks))


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch_geometry(output):
    file, metadata = output / "DIIID_geom.json", output / "geometry_provenance.json"
    if file.exists() and metadata.exists():
        provenance = json.loads(metadata.read_text(encoding="utf-8"))
        if sha256(file) != provenance["sha256"]:
            raise ValueError("cached geometry hash changed")
        return json.loads(file.read_text()), provenance
    req = urllib.request.Request("https://api.github.com/repos/hansec/OpenFUSIONToolkit/commits/main",
                                 headers={"User-Agent": "DIII-D-frozen-heat-research"})
    revision = json.load(urllib.request.urlopen(req, timeout=30))["sha"]
    url = f"https://raw.githubusercontent.com/hansec/OpenFUSIONToolkit/{revision}/{GEOMETRY_PATH}"
    content = urllib.request.urlopen(url, timeout=30).read()
    file.write_bytes(content)
    provenance = {"commit": revision, "url": url, "sha256": sha256(file),
                  "campaign_match": "unknown; public research baseline"}
    save_json(metadata, provenance)
    return json.loads(content), provenance


def load_case(shot, t_ms, method="cubic"):
    import pyarrow.parquet as pq
    columns = ["efit_times", "efit_psirz", "efit_grid_R", "efit_grid_Z",
               "efit_lcfs_n", "efit_lcfs_r", "efit_lcfs_z", "efit_r_axis", "efit_z_axis",
               "magnetics_bcoil", "magnetics_time"]
    s = pq.read_table(shot, columns=columns).to_pylist()[0]
    times = np.asarray(s["efit_times"])
    matches = np.flatnonzero(np.isclose(times, t_ms, atol=1e-6, rtol=0))
    if len(matches) != 1:
        raise ValueError("requested fixed time must exist exactly once")
    i = int(matches[0])
    tf = np.interp(t_ms, s["magnetics_time"], s["magnetics_bcoil"])
    field = FrozenField(s["efit_grid_R"], s["efit_grid_Z"], s["efit_psirz"][i], tf, method)
    n = int(s["efit_lcfs_n"][i])
    lcfs = np.column_stack([s["efit_lcfs_r"][i][:n], s["efit_lcfs_z"][i][:n]])
    axis = np.array([s["efit_r_axis"][i], s["efit_z_axis"][i]])
    return field, lcfs, axis


def separatrix(field, lcfs, axis):
    guess = lcfs[np.argmin(lcfs[:, 1])]
    solution = root(lambda p: field.gradient([p])[0], guess)
    xp = solution.x
    if (not solution.success or np.linalg.norm(xp-guess) > 0.15 or
            np.linalg.norm(field.gradient([xp])) > 1e-7):
        raise ValueError("lower X-point cannot be validated near LCFS")
    hessian = np.array([[field.spline.ev(xp[1], xp[0], dy=2),
                         field.spline.ev(xp[1], xp[0], dx=1, dy=1)],
                        [field.spline.ev(xp[1], xp[0], dx=1, dy=1),
                         field.spline.ev(xp[1], xp[0], dx=2)]])
    if np.linalg.det(hessian) >= 0:
        raise ValueError("stationary point is not a magnetic saddle")
    level = float(field.values([xp])[0])
    grid = np.linspace(axis[0]+0.03, field.r[-1]-0.1, 300)
    v = field.values(np.column_stack([grid, np.full(len(grid), axis[1])]))-level
    roots = [brentq(lambda r: field.values([[r, axis[1]]])[0]-level, a, b)
             for a, b, va, vb in zip(grid[:-1], grid[1:], v[:-1], v[1:]) if va*vb < 0]
    if len(roots) != 1:
        raise ValueError("unambiguous outboard separatrix crossing required")
    return roots[0], xp, level


def pretrace(field, r_sep, z_up, wall, tubes, step):
    edges = 0.1*np.linspace(0, 1, tubes+1)**2
    centers = (edges[:-1]+edges[1:])/2
    starts = np.column_stack([r_sep+centers, np.full(tubes, z_up)])
    paths = trace_paths(field, starts, step=step)
    # Starting outside a material wall is invalid input, even if it later enters.
    inside = shapely.contains(Polygon(wall), shapely.points(starts))
    lines = [LineString(p) for p in paths]
    hits = first_hits(lines, wall)
    hits.segment[~inside] = -1
    hits.distance[~inside] = np.inf
    hits.points[~inside] = np.nan
    drift = max(float(np.max(np.abs(field.values(p)-field.values([p[0]])[0]))) for p in paths)
    return edges, lines, inside, hits, drift


class ShapeEvaluator:
    def __init__(self, field, lcfs, wall, paths, inside, edges, r_sep, z_up, target_index, target_nodes=161):
        self.field, self.lcfs, self.wall = field, lcfs, wall
        self.paths, self.inside = paths, inside
        self.sources = [exponential_source(field, r_sep, z_up, edges, mm/1000) for mm in WIDTHS_MM]
        self.target_index = target_index
        self.target = np.linspace(wall[target_index], wall[target_index+1], target_nodes)
        self.target_length = np.linalg.norm(wall[target_index+1]-wall[target_index])
        # Cache all fixed segments without introducing a closing chord across target.
        fixed_paths = [LineString(wall[:target_index+1]), LineString(wall[target_index+1:])]
        self.fixed_geometry = shapely.union_all(fixed_paths)
        self.fixed_hits = self._fixed_hits()
        self.baseline = self.evaluate(np.zeros(4), baseline=True)
        self.baseline_peaks = np.array([max(m.q) for m in self.baseline["maps"]])
        self.evaluations = []
        self.best = None

    def _fixed_hits(self):
        # A temporary disjoint wall geometry is intersected; component indices
        # are identified against the original segments, excluding the target.
        from .frozen_mapping import _intersection_points
        intersections = shapely.intersection(np.asarray(self.paths, dtype=object), self.fixed_geometry)
        ids = np.full(len(self.paths), -1, dtype=int)
        pts = np.full((len(self.paths), 2), np.nan)
        ds = np.full(len(self.paths), np.inf)
        segments = np.diff(self.wall, axis=0)
        valid_segments = np.arange(len(segments)) != self.target_index
        for i, (line, inter) in enumerate(zip(self.paths, intersections)):
            if not self.inside[i]:
                continue
            candidates = _intersection_points(inter)
            if not candidates:
                continue
            lengths = [line.project(shapely.Point(p)) for p in candidates]
            k = int(np.argmin(lengths))
            point = candidates[k]
            a = np.clip(np.sum((point-self.wall[:-1])*segments, axis=1)/np.sum(segments**2, axis=1), 0, 1)
            residual = np.linalg.norm(self.wall[:-1]+a[:, None]*segments-point, axis=1)
            residual[~valid_segments] = np.inf
            ids[i], pts[i], ds[i] = np.argmin(residual), point, lengths[k]
        return Hits(ids, pts, ds)

    def evaluate(self, coefficients, baseline=False):
        target = deform_target(self.target, coefficients)
        length = float(np.sum(np.linalg.norm(np.diff(target, axis=0), axis=1)))
        if (length > 1.2*self.target_length or not LineString(target).is_simple or
                LineString(target[1:-1]).intersects(self.fixed_geometry) or
                LineString(target).intersects(Polygon(self.lcfs))):
            return {"geometry_valid": False, "feasible": False, "penalty": 1000.0,
                    "target": target, "length": length}
        candidate_wall = np.vstack([self.wall[:self.target_index], target, self.wall[self.target_index+2:]])
        target_segments = np.arange(self.target_index, self.target_index+len(target)-1)
        target_hits = first_hits(self.paths, target)
        choice = self.inside & (target_hits.distance < self.fixed_hits.distance)
        segment = self.fixed_hits.segment.copy()
        # fixed original segment indices after target shift by inserted segment count
        segment[segment > self.target_index] += len(target)-2
        segment[choice] = self.target_index+target_hits.segment[choice]
        points = self.fixed_hits.points.copy()
        points[choice] = target_hits.points[choice]
        hits = Hits(segment, points, np.minimum(target_hits.distance, self.fixed_hits.distance))
        maps = [map_heat(self.field, hits, s.weights, s.h, candidate_wall) for s in self.sources]
        # A tube with positive source power in ANY prescribed width is active.
        # Never remove a powered grazing hit by thresholding its surface q.
        active_power = np.max([m.hit_power for m in maps], axis=0)
        angle_status = maps[0].angle_application(ANGLE_FLOOR, target_segments, active_power)
        angle = angle_status["minimum_angle_deg"]
        mapped_to_target = [float(np.sum(m.segment_power[target_segments])) for m in maps]
        missing = max(m.unmapped_power for m in maps)
        power_error = max(m.balance_error for m in maps)
        feasible = candidate_feasible(angle, missing, power_error)
        # These are constraint violations, not changes to the prescribed source.
        penalties = constraint_penalties(angle_status, missing)
        penalty = sum(penalties.values())
        return {"geometry_valid": True, "feasible": feasible, "penalty": penalty,
                "target": target, "wall": candidate_wall, "maps": maps,
                "target_segments": target_segments, "angle": angle, "length": length,
                "target_power": mapped_to_target, "missing": missing,
                "angle_application": angle_status, "constraint_penalties": penalties}

    def objective(self, coefficients):
        result = self.evaluate(coefficients)
        if not result["geometry_valid"]:
            return result["penalty"]
        ratios = np.array([max(m.q) for m in result["maps"]])/self.baseline_peaks
        score = float(max(ratios))
        self.evaluations.append({"coefficients_m": list(coefficients), "worst_peak_ratio": score,
                                 "min_angle_deg": result["angle"], "unmapped_power": result["missing"],
                                 "angle_status": result["angle_application"]["status"],
                                 "constraint_penalties": result["constraint_penalties"],
                                 "feasible": result["feasible"]})
        if result["feasible"] and (self.best is None or score < self.best[0]):
            self.best = (score, np.asarray(coefficients).copy(), result)
        return score+result["penalty"]


def rows_for(evaluator, result, label):
    rows = []
    for mm, base_peak, source, m, target_power in zip(WIDTHS_MM, evaluator.baseline_peaks,
                                                   evaluator.sources, result["maps"], result["target_power"]):
        selected = np.isin(m.segment, result["target_segments"])
        rows.append({"geometry": label, "lambda_q_mm": int(mm), "q_over_P_peak_m2inv": float(max(m.q)),
                     "all_wall_peak_ratio": float(max(m.q)/base_peak),
                     "target_peak_m2inv": float(max(m.q[selected])) if np.any(selected) else 0,
                     "minimum_target_angle_deg": result["angle"], "target_power": target_power,
                     "angle_status": result["angle_application"]["status"],
                     "other_wall_power": m.mapped_power-target_power, "unmapped_power": m.unmapped_power,
                     "power_balance_error": m.balance_error, "feasible": result["feasible"]})
    return rows


def angle_scope_rows(evaluator, result, label):
    """Per-tube audit of the exact mask used for the candidate angle constraint."""
    heat = result["maps"][0]
    active_power = np.max([m.hit_power for m in result["maps"]], axis=0)
    applied = heat.angle_mask(result["target_segments"], active_power)
    target_hit = np.isin(heat.segment, result["target_segments"])
    rows = []
    for i, segment in enumerate(heat.segment):
        original = int(segment)
        if target_hit[i]:
            original = evaluator.target_index
        elif segment > result["target_segments"][-1]:
            original -= len(result["target"])-2
        reason = ("" if applied[i] else "no_first_hit" if segment < 0 else
                  "non_target_surface" if not target_hit[i] else "zero_source_power")
        row = {"geometry": label, "tube": i, "source_x_mm": float(evaluator.sources[0].x[i]*1000),
               "original_wall_segment": original, "is_explicit_target": bool(target_hit[i]),
               "first_hit_R_m": float(heat.points[i, 0]) if segment >= 0 else None,
               "first_hit_Z_m": float(heat.points[i, 1]) if segment >= 0 else None,
               "angle_deg": float(heat.angles[i]) if np.isfinite(heat.angles[i]) else None,
               "angle_floor_applied": bool(applied[i]), "excluded_reason": reason,
               "angle_lower_bound_violation_deg": float(max(0, ANGLE_FLOOR-heat.angles[i])) if applied[i] and np.isfinite(heat.angles[i]) else None}
        for mm, m in zip(WIDTHS_MM, result["maps"]):
            row[f"first_hit_power_lambda{mm}mm"] = float(m.hit_power[i])
        rows.append(row)
    return rows


def write_angle_scope(output, evaluator, result):
    rows = angle_scope_rows(evaluator, evaluator.baseline, "baseline")+angle_scope_rows(evaluator, result, "candidate")
    with (output/"angle_scope.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def plots(output, field, lcfs, xp, evaluator, result):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 7), constrained_layout=True)
    ax = axes[0]
    ax.contour(field.r, field.z, field.psi, levels=35, colors="0.85", linewidths=0.5)
    ax.plot(*np.vstack([lcfs, lcfs[0]]).T, color="tab:blue", label="EFIT LCFS")
    ax.plot(*evaluator.wall.T, color="0.25", label="Public PFC contour")
    ax.plot(*evaluator.target.T, color="tab:orange", linewidth=3, label="Baseline target")
    ax.plot(*result["target"].T, color="tab:red", linewidth=2, label="Candidate target")
    ax.plot(*xp, "x", color="tab:blue")
    for path in evaluator.paths[::max(1, len(evaluator.paths)//18)]:
        ax.plot(*np.array(path.coords).T, color="tab:green", alpha=0.35, linewidth=0.6)
    ax.set(xlabel="R (m)", ylabel="Z (m)", title="Frozen EFIT / first-hit target window")
    ax.set_aspect("equal"); ax.legend(fontsize=8)
    ax = axes[1]
    ax.plot(*evaluator.wall.T, color="0.6")
    ax.plot(*evaluator.target.T, color="tab:orange", label="Baseline")
    ax.plot(*result["target"].T, color="tab:red", label="Candidate")
    lo, hi = np.min(evaluator.target, axis=0), np.max(evaluator.target, axis=0)
    ax.set(xlim=(lo[0]-0.05, hi[0]+0.05), ylim=(lo[1]-0.08, hi[1]+0.08),
           xlabel="R (m)", ylabel="Z (m)", title="Target window; minimum full-field angle = 1.5 deg")
    ax.set_aspect("equal"); ax.legend()
    fig.savefig(output/"geometry.png", dpi=180); plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, mm, baseline, candidate in zip(axes.flat, WIDTHS_MM, evaluator.baseline["maps"], result["maps"]):
        for m, res, label in [(baseline, evaluator.baseline, "Baseline"), (candidate, result, "Candidate")]:
            selected = np.isin(m.segment, res["target_segments"])
            p, q = m.points[selected], m.q[selected]
            order = np.argsort(p[:, 0])
            ax.plot(p[order, 0], q[order], label=label)
        ax.set(title=f"lambda_q = {mm} mm", xlabel="R at target hit (m)", ylabel="q_surface / P_branch (1/m²)")
        ax.legend(fontsize=8)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0, 0.8, "Unit prescribed outer-branch power\nFixed equilibrium / direct mapping\nNo measured MW/m² claim\nAngle floor: 1.5 deg only", fontsize=11)
    fig.savefig(output/"heat_profiles.png", dpi=180); plt.close(fig)


def verify_candidate(shot, t_ms, wall, target_index, coefficients, tubes, step):
    """Hold angle=1.5 and design fixed while checking numerical representation."""
    checks = []
    for name, method, count, ds, nodes in [
            ("cubic_coarse", "cubic", tubes, step, 161),
            ("cubic_refined", "cubic", 2*tubes, step/2, 321),
            ("quintic_refined", "quintic", 2*tubes, step/2, 321)]:
        field, lcfs, axis = load_case(shot, t_ms, method)
        r_sep, _, _ = separatrix(field, lcfs, axis)
        edges, paths, inside, hits, drift = pretrace(field, r_sep, axis[1], wall, count, ds)
        evaluator = ShapeEvaluator(field, lcfs, wall, paths, inside, edges, r_sep, axis[1], target_index, nodes)
        result = evaluator.evaluate(coefficients)
        if not result["geometry_valid"]:
            checks.append({"name": name, "geometry_valid": False})
            continue
        ratios = [float(max(m.q)/b) for m, b in zip(result["maps"], evaluator.baseline_peaks)]
        checks.append({"name": name, "method": method, "tubes": count, "step_m": ds,
                       "surface_nodes": nodes, "angle_floor_deg": ANGLE_FLOOR,
                       "minimum_target_angle_deg": result["angle"], "feasible": result["feasible"],
                       "angle_application": result["angle_application"],
                       "constraint_penalties": result["constraint_penalties"],
                       "baseline_peaks": evaluator.baseline_peaks.tolist(),
                       "candidate_peaks": [float(max(m.q)) for m in result["maps"]],
                       "peak_ratios": ratios, "unmapped_power": result["missing"],
                       "max_path_psi_drift": drift,
                       "power_balance_error": max(m.balance_error for m in result["maps"])})
    return checks


def write_report(output, manifest, checks):
    candidates = [r for r in manifest["results"] if r["geometry"] == "candidate"]
    baseline = [r for r in manifest["results"] if r["geometry"] == "baseline"]
    feasible = verified_feasibility(manifest["feasible_candidate"], checks)
    if not checks:
        conclusion = "기준 형상 사전 확인; 후보 최적화·세분화 검증 전" if manifest["preflight_only"] else "저장 후보의 세분화 검증 전"
    elif manifest["angle_application"]["status"] == "target_unilluminated":
        conclusion = "지정 타깃에 첫 충돌점 없음; 입사각 위반과 구분"
    else:
        conclusion = "고정 제약을 만족한 연구 후보" if feasible else "저장 후보가 고정 제약의 수치 검증을 통과하지 못함"
    lines = ["# DIII-D 고정 평형 열유속 형상 탐색", "",
             f"입사각 하한은 **1.5° 하나**이며 입사각 민감도 계산은 수행하지 않았습니다.", "",
             f"선택한 좌표 구간에 대한 수치 결과: **{conclusion}**.", "",
             "공개 limiter 좌표에는 부품 태그가 없으므로 외부 설정으로 타깃 구간을 명시하고 geometry hash와 끝점 좌표를 확인합니다. 실제 샷의 DIII-D outer target 대응이 검증되지 않았다면 공개 형상의 연구 기준으로만 해석합니다.", "",
             f"고정 EFIT 시점은 {manifest['time_ms']:g} ms입니다. 공개 PFC 형상의 명시적으로 지정한 구간에 끝점과 끝점 기울기를 고정한 4변수 spline을 적용했습니다. 허용 변위 ±20 mm, 길이 증가는 20% 이내입니다.", "",
             f"기준 형상의 최소 전체 자기장 입사각: {angle_text(baseline[0]['minimum_target_angle_deg'], baseline[0].get('angle_status'))}; 저장한 후보: {angle_text(candidates[0]['minimum_target_angle_deg'], candidates[0].get('angle_status'))}. 타깃에 양의 파워를 가진 첫 충돌점이 없으면 각도를 정의하지 않으며 0° 위반으로 보고하지 않습니다.", "",
             "허용되지 않은 후보는 최적 형상으로 채택하지 않습니다. 탐색 실패가 이 평형에서 모든 가능한 디버터 설계가 불가능하다는 증명은 아닙니다. 현재 결과에는 표면 독립 적분으로 검증한 파워 잔차와 illuminated-support 상위 꼬리·젖음 면적 계산이 포함되지 않아 검증된 열유속 최적값으로 인증하지 않습니다.", "",
             "| λq (mm) | 기준 peak q/P (m⁻²) | 후보 peak q/P (m⁻²) | 후보/기준 | 타깃 흡수 파워 | 다른 벽 파워 |", "|---:|---:|---:|---:|---:|---:|"]
    for b, c in zip(baseline, candidates):
        lines.append(f"| {c['lambda_q_mm']} | {b['q_over_P_peak_m2inv']:.4f} | {c['q_over_P_peak_m2inv']:.4f} | {c['all_wall_peak_ratio']:.4f} | {c['target_power']:.6f} | {c['other_wall_power']:.6f} |")
    lines += ["", "## 모델과 단위", "",
              "파워는 **outer SOL branch 유입 파워 1**입니다. 전체 장치 P_SOL이나 검증된 inner/outer power split을 의미하지 않습니다. 표면 열유속 q/P 단위는 m⁻²입니다. 실제 MW/m²를 보고하지 않습니다.", "",
              "상류 q_parallel=C exp(-x/λq)를 공통 0–100 mm 범위에서 한 번 정규화하고, h=q_parallel/B_tot 및 q_surface=h|B·n|로 첫 표면에만 열을 배정했습니다. 미포착 열을 제거한 뒤 재정규화하지 않았습니다. 다른 벽의 열도 목적함수와 파워 수지에 포함됩니다.", "",
              f"최대 미포착 파워는 {max(r['unmapped_power'] for r in candidates):.3g}, 전체 파워 분할 수지 잔차는 {max(r['power_balance_error'] for r in candidates):.3g}입니다. 이 분할 수지는 표면 국소 열유속 적분의 독립적인 실측 검증을 뜻하지 않습니다.", "",
              "전체 입사각은 asin(|B·n|/|B|)입니다. B_phi는 공개 총 144권선과 bcoil 전류로 계산한 이상 축대칭 진공 근사이며 토로이달 ripple과 calibration 불확실성을 포함하지 않습니다.", "",
              "1.5°는 DIII-D 문헌의 대표 하부 shelf 입사각을 사용자가 연구 하한으로 채택한 값입니다. 공식 보편적 제작·안전 한계라고 주장하지 않습니다.", "",
              "## 수치 확인", "",
              "| 표현 | tubes | 최소 입사각 (°) | 최대 peak ratio | 미포착 파워 | 최대 ψ 경로 drift | 허용 |", "|---|---:|---:|---:|---:|---:|---|"]
    for c in checks:
        if "peak_ratios" in c:
            lines.append(f"| {c['name']} | {c['tubes']} | {angle_text(c['minimum_target_angle_deg'], c.get('angle_application', {}).get('status'))} | {max(c['peak_ratios']):.4f} | {c['unmapped_power']:.3g} | {c['max_path_psi_drift']:.3g} | {c['feasible']} |")
    lines += ["", "원본 ψ 격자는 약 ΔR=26.6 mm, ΔZ=50 mm입니다. mm 규모 λq의 계산은 보간된 평형을 가정한 연구이며 실제 장치의 mm 수준 정확도를 검증하지 않습니다. 비교에는 같은 원본 ψ의 cubic 및 quintic 보간과 각 보간의 해석적 미분을 사용하며, 경로 ψ drift도 함께 보고합니다.", "",
              "geometry.json은 샷의 실제 캠페인과 일치한다고 확인되지 않았습니다. 재료, 방사, 중성입자, sheath 열전달, cross-field spreading, detachment 및 ELM은 이 직접 매핑 모델에 포함하지 않습니다.", "",
              "## 재현 및 자료", "", "입력 hash·geometry commit·설정은 `manifest.json`, 전체 평가 이력은 `evaluations.json`, 수치 비교는 `numerical_checks.json`에 있습니다.", "",
              "- [공개 geometry 예제](https://openfusiontoolkit.github.io/OpenFUSIONToolkit/docs/v26.6/doc_tMaker_DIIID_ex1.html)",
              "- [DIII-D 입사각·HEAT 연구](https://www.tandfonline.com/doi/full/10.1080/15361055.2025.2493406)",
              "- [TF 권선 문헌](https://fusion.gat.com/pubs-ext/SOFT02/A24059.pdf)",
              "- [열유속 폭 모델](https://courses.physics.ucsd.edu/2021/Spring/physics218c/AA_Eich_2013_Nucl._Fusion_53_093031.pdf)", ""]
    if manifest.get("target_configuration"):
        config = manifest["target_configuration"]
        lines += ["## 명시적 타깃 설정", "",
                  f"설정 이름: {config['name']}; 원본 segment: {config['segment_index']}; 끝점 R,Z (m): {config['endpoints_RZ_m']}. 실제 하드웨어 정체 확인: {config['physical_identity_validated']}.", "",
                  "1.5°는 지정 타깃의 첫 충돌점에만 적용합니다. 다른 벽에 먼저 닿는 열은 그대로 다른 벽 파워와 전체 벽 peak 목적함수에 포함하며 타깃으로 옮기지 않습니다.", "",
                  f"타깃 판정 상태: {manifest['angle_application']['status']}; 첫 충돌점 수: {manifest['angle_application']['target_hit_count']}. 타깃 미조사는 유효한 열유속 최적화 결과로 채택하지 않습니다.", ""]
        if "applied_hit_count" in manifest["angle_application"]:
            lines += [f"각도 하한 적용 점 수: {manifest['angle_application']['applied_hit_count']}; 소스 파워 0으로 제외된 타깃 교차점: {manifest['angle_application']['excluded_zero_power_target_hit_count']}.", "",
                      "적용 조건은 ‘명시적 타깃의 첫 충돌점이며, 폭 스캔 중 하나 이상의 유입 tube 파워가 양수’입니다. 임의의 열유속 cutoff로 낮은 각도를 제외하지 않습니다. 차폐되어 먼저 다른 벽에 닿는 경로는 비타깃으로 처리합니다.", "",
                      "점별 적용/제외 사유와 원본 벽 segment, 좌표 및 폭별 파워는 `angle_scope.csv`에 기록합니다. 타깃 미조사 벌점과 각도 하한 위반 벌점은 `constraint_penalties`에서 별도 항목으로 관리합니다.", ""]
    if manifest.get("angle_application_diagnostic", {}).get("minimum_angle_hit_RZ_m") is not None:
        d = manifest["angle_application_diagnostic"]
        lines += ["## 입사각 적용 위치 확인", "",
                  f"기준 형상에서 최소각 위치는 R={d['minimum_angle_hit_RZ_m'][0]:.6f} m, Z={d['minimum_angle_hit_RZ_m'][1]:.6f} m입니다. 해당 tube의 상류 separatrix 거리 x={d['minimum_angle_source_x_mm']:.6f} mm; λq=1 mm 전체 벽 최대 열유속 위치와 일치 여부: {d['minimum_angle_is_lambda1mm_peak_location']}.", "",
                  f"같은 평형의 R=1.47 m, Z=-1.2448 m 수평 shelf 점에서 전체 입사각은 {d['comparison_shelf_angle_deg']:.4f}°입니다. 이 점은 현재 선택된 열 침착 구간과 다르며, 해당 점에서 1.5°를 만족한다고 선택 구간의 최소각도 만족한다고 결론 내리지 않습니다.", "",
                  "각도 조건의 물리적 적용 대상을 확정하려면 지정한 좌표 구간과 실제 target/baffle/pump 구조 태그를 대응시켜야 합니다.", ""]
    (output/"REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shot", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-config", required=True, type=Path,
                        help="explicit target name, geometry hash, segment and endpoint coordinates")
    parser.add_argument("--time-ms", type=float, default=1400)
    parser.add_argument("--tubes", type=int, default=320)
    parser.add_argument("--iterations", type=int, default=25)
    parser.add_argument("--step", type=float, default=0.002)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true", help="re-evaluate the saved fixed design; no angle scan")
    args = parser.parse_args(argv)
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    geom, provenance = fetch_geometry(output)
    wall = np.asarray(geom["limiter"], dtype=float)
    wall = np.vstack([wall, wall[0]])
    target_config = json.loads(args.target_config.read_text(encoding="utf-8"))
    target_index = select_target(wall, provenance["sha256"], target_config)
    if args.verify_only:
        manifest = json.loads((output/"manifest.json").read_text(encoding="utf-8"))
        if manifest["shot_sha256"] != sha256(args.shot):
            raise ValueError("verification shot must match the original case hash")
        if manifest.get("target_configuration") != target_config:
            raise ValueError("verification requires the same explicit target configuration; rerun legacy automatic cases")
        if manifest.get("angle_constraint_scope") != ANGLE_SCOPE:
            raise ValueError("angle scope policy changed; rerun the case before verification")
        checks = verify_candidate(args.shot, manifest["time_ms"], wall,
                                  manifest["target_original_segment"], manifest["coefficients_m"],
                                  manifest["tubes"], manifest["step_m"])
        save_json(output/"numerical_checks.json", checks)
        field, lcfs, axis = load_case(args.shot, manifest["time_ms"])
        r_sep, _, _ = separatrix(field, lcfs, axis)
        edges, paths, inside, _, _ = pretrace(field, r_sep, axis[1], wall, manifest["tubes"], manifest["step_m"])
        e = ShapeEvaluator(field, lcfs, wall, paths, inside, edges, r_sep, axis[1], manifest["target_original_segment"])
        result = e.evaluate(manifest["coefficients_m"])
        if not result["geometry_valid"]:
            raise ValueError("stored candidate geometry is invalid; cannot audit its angle constraint")
        manifest["angle_application"] = result["angle_application"]
        manifest["constraint_penalties"] = result["constraint_penalties"]
        write_angle_scope(output, e, result)
        heat = e.baseline["maps"][0]
        selected = heat.angle_mask(e.baseline["target_segments"], np.max([m.hit_power for m in e.baseline["maps"]], axis=0))
        minimum = int(np.flatnonzero(selected)[np.argmin(heat.angles[selected])]) if np.any(selected) else None
        manifest["angle_application_diagnostic"] = {
            "application": "first-hit points on explicitly configured target only",
            "physical_target_identity_validated": target_config["physical_identity_validated"],
            "minimum_angle_hit_RZ_m": heat.points[minimum].tolist() if minimum is not None else None,
            "minimum_angle_source_x_mm": float(e.sources[0].x[minimum]*1000) if minimum is not None else None,
            "minimum_angle_is_lambda1mm_peak_location": bool(minimum == int(np.argmax(heat.q))) if minimum is not None else None,
            "comparison_shelf_RZ_m": [1.47, -1.2448],
            "comparison_shelf_angle_deg": float(field.incidence([[1.47,-1.2448]], [[0,1]])[0])}
        manifest["numerically_verified_feasible"] = verified_feasibility(manifest["feasible_candidate"], checks)
        manifest["certified_heat_optimum"] = False
        manifest["surface_integral_power_residual"] = None
        save_json(output/"manifest.json", manifest)
        write_report(output, manifest, checks)
        print(json.dumps(checks, indent=2), flush=True)
        return
    field, lcfs, axis = load_case(args.shot, args.time_ms)
    r_sep, xp, psi_sep = separatrix(field, lcfs, axis)
    print(f"Fixed time {args.time_ms:g} ms; OMP R_sep={r_sep:.6f}; lower X-point={xp}", flush=True)
    edges, paths, inside, hits, drift = pretrace(field, r_sep, axis[1], wall, args.tubes, args.step)
    sources = [exponential_source(field, r_sep, axis[1], edges, mm/1000) for mm in WIDTHS_MM]
    print(f"Explicit target {target_config['name']} / segment {target_index}: {wall[target_index:target_index+2].tolist()}; flux drift={drift:.3e}", flush=True)
    evaluator = ShapeEvaluator(field, lcfs, wall, paths, inside, edges, r_sep, axis[1], target_index)
    if not evaluator.baseline["geometry_valid"]:
        raise ValueError("public baseline target intersects EFIT LCFS or fixed wall; cannot optimize this window")
    print(f"Baseline min angle={angle_text(evaluator.baseline['angle'], evaluator.baseline['angle_application']['status'])}; status={evaluator.baseline['angle_application']['status']}; max unmapped={evaluator.baseline['missing']:.5g}", flush=True)
    coefficients = np.zeros(4)
    optimizer_result = None
    if not args.preflight_only:
        # Track the zero design as a candidate too; never worsen a feasible baseline.
        evaluator.objective(coefficients)
        optimum = differential_evolution(evaluator.objective, [(-0.02, 0.02)]*4,
                                         seed=42, maxiter=args.iterations, popsize=7,
                                         polish=False, tol=0.003, workers=1)
        optimizer_result = {"success": bool(optimum.success), "message": str(optimum.message),
                            "nfev": int(optimum.nfev), "objective_with_penalty": float(optimum.fun)}
        coefficients = evaluator.best[1] if evaluator.best is not None else optimum.x
    result = evaluator.evaluate(coefficients)
    if not result["geometry_valid"]:
        result, coefficients = evaluator.baseline, np.zeros(4)
    rows = rows_for(evaluator, evaluator.baseline, "baseline")+rows_for(evaluator, result, "candidate")
    with (output/"width_scan.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    np.savez_compressed(output/"case_arrays.npz", R=field.r, Z=field.z, psi=field.psi, lcfs=lcfs,
                        baseline_wall=wall, candidate_wall=result["wall"], target=result["target"],
                        coefficients=coefficients, x_sources=sources[0].x)
    manifest = {"shot_path": str(args.shot.resolve()), "shot_sha256": sha256(args.shot),
                "geometry": provenance, "time_ms": args.time_ms, "angle_floor_deg": ANGLE_FLOOR,
                "angle_sensitivity": False, "widths_mm": WIDTHS_MM.tolist(), "outer_branch_power": 1.0,
                "source_max_x_m": 0.1, "source_profile": "q_parallel,u proportional to exp(-x/lambda_q)",
                "tf_model": "ideal axisymmetric 144-turn field, bcoil kA to A", "tf_current_ka": field.tf_current_ka,
                "interpolation": field.method, "r_sep_m": r_sep, "xpoint_m": xp.tolist(),
                "psi_sep_Wb_per_rad": psi_sep, "max_path_psi_drift": drift,
                "target_original_segment": target_index, "tubes": args.tubes, "step_m": args.step,
                "target_configuration": target_config, "angle_application": result["angle_application"],
                "angle_constraint_scope": ANGLE_SCOPE, "constraint_penalties": result["constraint_penalties"],
                "python": platform.python_version(), "elapsed_seconds": time.monotonic()-started,
                "preflight_only": args.preflight_only, "optimizer": optimizer_result,
                "feasible_candidate": result["feasible"], "coefficients_m": coefficients.tolist(),
                "minimum_target_angle_deg": result["angle"], "results": rows,
                "numerically_verified_feasible": False, "certified_heat_optimum": False,
                "scope": "conditional direct magnetic heat mapping; unknown campaign match; no measured heat load"}
    save_json(output/"manifest.json", manifest)
    save_json(output/"evaluations.json", evaluator.evaluations)
    save_json(output/"numerical_checks.json", [])
    write_angle_scope(output, evaluator, result)
    plots(output, field, lcfs, xp, evaluator, result)
    write_report(output, manifest, [])
    status = case_status(result)
    print(f"{status}; output={output}; elapsed={manifest['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
