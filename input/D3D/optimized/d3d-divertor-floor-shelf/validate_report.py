"""Verify corrected components, report heat reduction, and show both targets."""
import csv
import json
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import substring
from divertor_optimization.d3d_multisurface import MultiSurfaceEvaluator
from divertor_optimization.frozen_mapping import deform_target
from divertor_optimization.d3d_frozen_case import load_case,separatrix,pretrace,save_json,sha256,WIDTHS_MM

p=Path(__file__).resolve().parent
m=json.loads((p/"manifest.json").read_text(encoding="utf-8"))
assert sha256(m["shot_path"])==m["shot_sha256"]
assert sha256(p/'DIIID_geom.json')==m['config']['geometry_sha256']
geometry=json.loads((p/"DIIID_geom.json").read_text())
wall=np.asarray(geometry["limiter"],dtype=float)
wall=np.vstack([wall,wall[0]])
vessel=shapely.union_all([Polygon(v["pts"]) for v in geometry["vv"]])
evaluators=[]
for name,method,count,step,nodes in [
    ("cubic_coarse","cubic",320,.002,161),("cubic_refined","cubic",640,.001,321),
    ("quintic_refined","quintic",640,.001,321)]:
    f,lcfs,axis=load_case(m["shot_path"],m["time_ms"],method)
    r_sep,xp,_=separatrix(f,lcfs,axis)
    edges,paths,inside,_,drift=pretrace(f,r_sep,axis[1],wall,count,step)
    e=MultiSurfaceEvaluator(f,lcfs,wall,paths,inside,edges,r_sep,axis[1],m["config"]["surfaces"],nodes,vessel_solids=vessel)
    evaluators.append((name,e,drift))
    print("prepared",name,flush=True)
records=[json.loads(s) for s in (p/"evaluations.jsonl").read_text(encoding="utf-8").splitlines()]
assert len({tuple(r["coefficients_m"]) for r in records})==len(records)
vessel_audit=[]
for row in records:
    safe=False
    if row['geometry_valid']:
        safe=all(not LineString(deform_target(np.linspace(wall[s['segment_index']],wall[s['segment_index']+1],321),
                    row['coefficients_m'][4*j:4*j+4])).intersects(vessel)
                 for j,s in enumerate(m['config']['surfaces']))
    row['vessel_safe_geometry_valid']=bool(row['geometry_valid'] and safe)
    vessel_audit.append({'case':row['case'],'search_geometry_valid':row['geometry_valid'],
                         'vessel_safe_geometry_valid':row['vessel_safe_geometry_valid']})
with (p/'vessel_geometry_audit.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(vessel_audit[0]));w.writeheader();w.writerows(vessel_audit)
m.setdefault('search_counts_before_vessel_audit',m['counts'])
m['vessel_constraint_note']='Stored search log precedes vessel-solid filter; all geometrically valid logged candidates are audited at 321 nodes. Current run_search.py enforces this constraint during future searches, so its rerun trajectory may differ.'
m['counts']={'evaluated':len(records),'geometry_valid':sum(r['vessel_safe_geometry_valid'] for r in records),
             'feasible':sum(r['feasible'] and r['vessel_safe_geometry_valid'] for r in records)}
assert m['counts']['geometry_valid']>=1000
ranked=sorted((r for r in records if r['feasible'] and r['vessel_safe_geometry_valid']),key=lambda r:r['worst_all_wall_peak_ratio'])
attempts=[]; chosen=None
for rank,row in enumerate(ranked):
    checks=[]; results=[]
    for name,e,drift in evaluators:
        r=e.evaluate(row["coefficients_m"])
        if not r["geometry_valid"]:
            checks.append({"name":name,"feasible":False,"geometry_valid":False}); break
        ratios=[float(max(h.q)/b) for h,b in zip(r["maps"],e.baseline_peaks)]
        checks.append({"name":name,"feasible":r["feasible"],"geometry_valid":True,
                       "peak_ratios":ratios,"angles":r["angles"],"unmapped_power":r["missing"],
                       "surface_power":r["surface_power"],"power_balance_error":r["power_balance_error"],
                       "angle_penalty":r["angle_penalty"],"max_path_psi_drift":drift})
        checks[-1]['vessel_clearance_m']={n:LineString(t).distance(vessel) for n,t in r['targets'].items()}
        results.append(r)
        if not r["feasible"]: break
    passed=len(checks)==3 and all(c["feasible"] for c in checks)
    attempts.append({"case":row["case"],"rank":rank,"passed":passed,"checks":checks})
    if rank%32==0: print("verification",rank,"passed",passed,flush=True)
    if passed:
        chosen=row,checks,results[0]; break
assert chosen is not None
row,checks,r=chosen; e=evaluators[0][1]
table=[]
for label,result in [("baseline",e.baseline),("candidate",r)]:
    for j,(mm,heat,bpeak) in enumerate(zip(WIDTHS_MM,result["maps"],e.baseline_peaks)):
        idx=int(np.argmax(heat.q))
        entry={"geometry":label,"lambda_q_mm":int(mm),"all_wall_peak_q_over_P":float(max(heat.q)),
               "all_wall_peak_ratio":float(max(heat.q)/bpeak),
               "peak_R_m":float(heat.points[idx,0]),"peak_Z_m":float(heat.points[idx,1]),
               "peak_original_segment":int(result["original_segment"][heat.segment[idx]]),
               "floor_peak_q_over_P":float(max(heat.q[np.isin(heat.segment,result["surface_segments"]["floor"])],default=0.)),
               "shelf_peak_q_over_P":float(max(heat.q[np.isin(heat.segment,result["surface_segments"]["shelf"])],default=0.)),
               "floor_power":result["surface_power"]["floor"][j],"shelf_power":result["surface_power"]["shelf"][j],
               "minimum_floor_angle_deg_diagnostic":result["angles"]["floor"]["minimum_angle_deg"],
               "minimum_shelf_angle_deg":result["angles"]["shelf"]["minimum_angle_deg"],
               "unmapped_power":heat.unmapped_power,"power_balance_error":heat.balance_error,"feasible":result["feasible"]}
        entry["fixed_surface_power"]=heat.mapped_power-entry["floor_power"]-entry["shelf_power"]
        table.append(entry)
with (p/"width_scan.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
audit=[]
for label,result in [("baseline",e.baseline),("candidate",r)]:
    heat=result["maps"][0]; design,angle=e.hit_masks(result)
    for i,segment in enumerate(heat.segment):
        name="fixed_surface"
        for n,ids in result["surface_segments"].items():
            if segment in ids: name=n
        if segment<0: name="unmapped"
        audit.append({"geometry":label,"tube":i,"first_hit_component":name,
                      "original_segment":int(result["original_segment"][segment]) if segment>=0 else -1,
                      "R_m":float(heat.points[i,0]) if segment>=0 else None,
                      "Z_m":float(heat.points[i,1]) if segment>=0 else None,
                      "is_design_hit":bool(design[i]),"angle_floor_applied":bool(angle[i]),
                      "angle_deg":float(heat.angles[i]) if np.isfinite(heat.angles[i]) else None,
                      "reason":("shelf_positive_power_first_hit" if angle[i] else
                                "floor_has_no_angle_floor" if name=="floor" else
                                "no_positive_power_shelf_hit" if name=="shelf" else name)})
assert all(a["first_hit_component"]=="shelf" for a in audit if a["angle_floor_applied"])
with (p/"angle_scope.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
np.savez_compressed(p/"case_arrays.npz",R=e.field.r,Z=e.field.z,psi=e.field.psi,lcfs=e.lcfs,
                    baseline_wall=wall,candidate_wall=r["wall"],baseline_floor=e.targets["floor"],
                    candidate_floor=r["targets"]["floor"],baseline_shelf=e.targets["shelf"],
                    candidate_shelf=r["targets"]["shelf"],coefficients=row["coefficients_m"])
m.update(coefficients_m=row["coefficients_m"],retained_case=row["case"],verification_attempt_count=len(attempts),
         candidate_angles=r["angles"],numerically_verified_feasible=True,results=table,
         worst_peak_ratio=row["worst_all_wall_peak_ratio"],peak_reduction_percent=100*(1-row["worst_all_wall_peak_ratio"]),
         checked_reduction_percent=[100*(1-max(c["peak_ratios"])) for c in checks],
         max_actual_displacement_m=row["max_displacement_m"])
save_json(p/"manifest.json",m); save_json(p/"numerical_checks.json",checks);save_json(p/"candidate_verification.json",attempts)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,2,figsize=(12,6.5),constrained_layout=True)
for ax in axes:
    for v in geometry['vv']:
        pts=np.asarray(v['pts']);ax.fill(*pts.T,color='.85',zorder=0)
    ax.plot(*wall.T,color=".65",label="Public reference contour")
    for n,color in [("floor","tab:red"),("shelf","tab:purple")]:
        ax.plot(*e.targets[n].T,"--",color=color,linewidth=2,label=f"Baseline {n}")
        ax.plot(*r["targets"][n].T,color=color,linewidth=2.5,label=f"Candidate {n}")
    ax.plot(*np.vstack([e.lcfs,e.lcfs[0]]).T,color="tab:blue",linewidth=1.4,label="Fixed EFIT LCFS")
    ax.set(xlabel="R (m)",ylabel="Z (m)");ax.set_aspect("equal")
ax=axes[0];ax.set(title="Fixed equilibrium / two design surfaces")
heat=r["maps"][2]
for i in range(0,len(e.paths),max(1,len(e.paths)//18)):
    if heat.segment[i]>=0:
        line=substring(e.paths[i],0,e.paths[i].project(Point(heat.points[i])))
        ax.plot(*np.array(line.coords).T,color="tab:green",alpha=.4,linewidth=.6)
ax=axes[1];ax.set(xlim=(1.,1.82),ylim=(-1.49,-1.10),title="Floor: no angle floor / Shelf: 1.5 deg only")
peak=table[7];ax.plot(peak["peak_R_m"],peak["peak_Z_m"],"*",color="black",markersize=10,label="Candidate peak at lambda=3 mm")
ax.legend(fontsize=8,loc="upper right")
fig.savefig(p/"geometry.png",dpi=180);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
for ax,j,mm in zip(axes.flat,range(5),WIDTHS_MM):
    for n,color in [("floor","tab:red"),("shelf","tab:purple")]:
        for result,style,label in [(e.baseline,"--","baseline"),(r,"-","candidate")]:
            h=result["maps"][j];selected=np.isin(h.segment,result["surface_segments"][n]);order=np.argsort(h.points[selected,0])
            ax.plot(h.points[selected,0][order],h.q[selected][order],style,color=color,label=f"{n} {label}")
    ax.set(title=f"lambda_q={mm} mm",xlabel="R at first hit (m)",ylabel="q / P_outer (1/m²)")
axes.flat[0].legend(fontsize=7);axes.flat[-1].axis("off")
axes.flat[-1].text(0,.7,"Both floor and shelf are design variables.\n1.5 deg floor applies ONLY to shelf.\nPower is unit outer-branch power.",fontsize=10)
fig.savefig(p/"heat_profiles.png",dpi=180);plt.close(fig)
counts=m["counts"]
lines=["# DIII-D 하부 디버터: 바닥 타깃 + shelf 최적화", "",
       "기존 shelf 전용 탐색에서 제외됐던 주 열부하 바닥 타깃을 설계 변수에 포함했습니다. 1.5°는 사용자 선택에 따라 shelf에만 적용하며 바닥 각도는 진단값입니다.", "",
       f"서로 다른 {counts['evaluated']:,}개 후보 중 기하학적으로 유효한 {counts['geometry_valid']:,}개를 평가했습니다. 기본 해상도에서 제약을 만족한 후보는 {counts['feasible']:,}개입니다. 최소 1,000개 유효 후보 조건을 충족했습니다.", "",
       f"유지한 후보의 전체 벽 worst-case peak 감소는 기본 계산에서 {m['peak_reduction_percent']:.3f}%입니다. 세분화·보간 확인 감소율은 {m['checked_reduction_percent']}%입니다. 원본 평가 case {row['case']}를 유지했습니다.", "",
       "## 설계 및 제약 범위", "",
       "공개 floor segment 68과 shelf segment 59를 동시에 변형합니다. 각각 4변수 spline으로 끝점·끝점 기울기를 고정하며 계수는 ±100 mm 범위, 표면 길이는 각각 기준 대비 120% 이하입니다. 자기교차, 고정 벽 침범, 두 설계 표면 교차, LCFS 침범 및 공개 진공용기 고체와의 교차는 금지합니다. 용기 검사는 전체 탐색 후 재분류했으며 vessel_geometry_audit.csv에 원래 판정과 수정 판정을 보존했습니다. 나머지 표면은 고정하되 전체 벽 열부하와 파워에 포함합니다. 비설계 표면이 실제 비디버터라는 의미는 아닙니다.", "",
       "평형은 원본 EFIT 1400 ms로 고정하고 폭은 1·2·3·5·10 mm, outer branch 유입 파워는 1입니다. 하한 민감도 스캔은 하지 않습니다. 각도 제약은 shelf의 양의 파워 첫 충돌점에만 적용하며, 바닥·차폐된 후속 교차점·가열되지 않는 점에는 적용하지 않습니다. shelf에 열이 닿지 않으면 각도 위반으로 처리하지 않으며 전체 설계 표면이 열을 받는지 별도로 확인합니다.", "",
       f"spline 계수 (m), 순서 floor 4개 + shelf 4개: {row['coefficients_m']}. 실제 최대 변위 (m): {row['max_displacement_m']}.", "",
       "## 전체 벽 열유속", "",
       "| λq (mm) | 기준 전체벽 peak q/P (m⁻²) | 후보 전체벽 peak q/P | 감소 (%) | 바닥 흡수 파워 | shelf 흡수 파워 | 고정 표면 파워 |",
       "|---:|---:|---:|---:|---:|---:|---:|"]
for b,c in zip(table[:5],table[5:]):
    lines.append(f"| {c['lambda_q_mm']} | {b['all_wall_peak_q_over_P']:.6g} | {c['all_wall_peak_q_over_P']:.6g} | {100*(1-c['all_wall_peak_ratio']):.3f} | {c['floor_power']:.6g} | {c['shelf_power']:.6g} | {c['fixed_surface_power']:.6g} |")
lines += ["", "## 동일 하한으로 수치 재평가", "",
          "| 표현 | shelf 최소각 (°) | 바닥 최소각: 진단만 (°) | 전체 벽 worst ratio | 미포착 파워 |", "|---|---:|---:|---:|---:|"]
for c in checks:
    shelf=c["angles"]["shelf"]["minimum_angle_deg"];floor=c["angles"]["floor"]["minimum_angle_deg"]
    lines.append(f"| {c['name']} | {shelf if shelf is not None else 'N/A'} | {floor if floor is not None else 'N/A'} | {max(c['peak_ratios']):.6f} | {c['unmapped_power']:.3g} |")
lines += ["", "## 해석 한계", "",
          "공개 윤곽의 floor/shelf 명칭은 공개 하부 디버터 문헌과 좌표·고정 평형의 열 연결을 대응한 연구 기준입니다. 샷 캠페인·부품 정체는 미확인이며 제작 가능한 형상으로 인증하지 않습니다. 원본 ψ 격자는 약 26.6×50 mm입니다. 실제 MW/m², sheath·방사·중성입자·cross-field spreading·detachment·ELM·재료 온도를 재현하지 않습니다. 파워를 첫 교차에만 배정하고 포착된 파워를 재정규화하지 않았습니다. 파워 분할 수지는 독립적인 표면 q 적분 검증을 대신하지 않습니다. 결과는 이 수치 모델에 조건부인 후보이며 인증된 열유속 최적값은 아닙니다.", "",
          "재현: run_search.py와 validate_report.py. 모든 후보 evaluations.jsonl; 설정과 hash manifest.json; 점별 설계/각도 마스크 angle_scope.csv; 수치 확인 numerical_checks.json."]
(p/"REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
import runpy
runpy.run_path(str(p/'render_geometry.py'),run_name='__main__')
print("VERIFIED",m["peak_reduction_percent"],m["checked_reduction_percent"],m["counts"],flush=True)
