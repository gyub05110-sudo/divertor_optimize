"""Static research figures and exact tables from the refined OMP study."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from d3d_optimization.broadening_study import save,load_json
from d3d_optimization.tangent_geometry import Curve

COLORS={20:'#0072B2',50:'#D55E00',100:'#009E73'}


def run(folder):
    study=load_json(folder/'study.json'); base=np.load(folder/'baseline_refined.npz')
    data={mm:np.load(folder/f'R{mm:03d}'/'selected_refined.npz') for mm in COLORS}
    vals={mm:load_json(folder/f'R{mm:03d}'/'validation.json') for mm in COLORS}
    manifests={mm:load_json(folder/f'R{mm:03d}'/'manifest.json') for mm in COLORS}
    if not all(v['complete'] for v in vals.values()):raise ValueError('validation must finish before visualization')
    rows=[]
    for mm,v in vals.items():
        level=v['levels'][2];b=level['rows']['baseline'];s=level['rows']['selected'];g=level['geometry_audit']
        m=manifests[mm]
        rows.append(dict(radius_constraint_mm=mm,baseline_peak_MW_m2=b['peak_MW_m2'],peak_MW_m2=s['peak_MW_m2'],
            reduction_percent=100*(1-s['peak_ratio']),actual_floor_minimum_radius_mm=1000*g['actual_minimum_radius_m']['floor'],
            minimum_gap_mm=1000*g['gap_m'],required_gap_mm=1000*g['gap_minimum_m'],tangent_error_deg=g['tangent_angle_error_deg'],
            resolution_screen_passed=v['verdict']['resolution_screen_passed'],improvement_confirmed=v['verdict']['improvement_confirmed'],
            power_deposited_W=s['ledger']['deposited']*study['reference']['power_W'],power_unhit_W=s['ledger']['unhit']*study['reference']['power_W'],
            power_balance_error=s['ledger']['balance_error'],length_ratio=s['precise_length_ratios']['floor'],
            L_parallel_mean_m=s['thermal']['length_m']['power_weighted_mean'],sigma_mean_mm=1000*s['sigma_physical_mean_m'],
            search_evaluations=m['search']['evaluated'],feasible_evaluations=m['search']['feasible_evaluations'],
            search_seconds=m['search']['elapsed_seconds']))
    save(folder/'summary.json',rows)
    with (folder/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    power=study['reference']['power_W']
    def cells(a):
        wall=a['receiving_wall'];q=a['cell_power']/np.diff(a['area_edges'])*power/1e6
        return wall,q
    maximum=max([rows[0]['baseline_peak_MW_m2']]+[r['peak_MW_m2'] for r in rows])
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for ax,key in zip(axes.flat,[None,20,50,100]):
        a=base if key is None else data[key];wall,q=cells(a)
        lc=LineCollection(np.stack([wall[:-1],wall[1:]],axis=1),cmap='inferno',norm=plt.Normalize(0,maximum),linewidth=4)
        lc.set_array(q);ax.add_collection(lc)
        peak=rows[0]['baseline_peak_MW_m2'] if key is None else next(r['peak_MW_m2'] for r in rows if r['radius_constraint_mm']==key)
        name='Original baseline' if key is None else f'G1 + Rmin {key} mm + opening constraint'
        ax.set(xlim=(1.07,1.48),ylim=(-1.44,-1.23),xlabel='R [m]',ylabel='Z [m]',title=f'{name}\nPeak {peak:.4f} MW/m²')
        ax.set_aspect('equal');ax.grid(alpha=.15)
    fig.colorbar(lc,ax=axes,label='Cell-average heat flux [MW/m²]',shrink=.85)
    fig.suptitle('OMP-only geometry search | 1 MW source | same baseline and transport inputs',fontsize=15)
    fig.savefig(folder/'heat_flux_comparison.png',dpi=180);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(16,5.5),layout='constrained')
    for ax in axes:
        ax.plot(base['wall'][:,0],base['wall'][:,1],'--',color='#555555',lw=1.6,label='Original baseline')
        for mm,a in data.items():ax.plot(a['wall'][:,0],a['wall'][:,1],color=COLORS[mm],lw=1.7,label=f'Rmin {mm} mm')
        ax.set_aspect('equal');ax.grid(alpha=.18);ax.set(xlabel='R [m]',ylabel='Z [m]')
    axes[0].set(xlim=(1.10,1.45),ylim=(-1.44,-1.27),title='OMP floor and fixed neighbours')
    axes[1].set(xlim=(1.135,1.215),ylim=(-1.39,-1.32),title='Left junction: same tangent, fixed point')
    axes[2].set(xlim=(1.35,1.405),ylim=(-1.395,-1.315),title='Right opening: tip-to-floor clearance')
    for mm,v in vals.items():
        g=v['levels'][2]['geometry_audit'];p=np.array([g['lip_point'],g['gap_floor_point']])
        axes[2].plot(p[:,0],p[:,1],':',color=COLORS[mm],lw=1.5)
    axes[0].legend(loc='lower left',fontsize=9)
    fig.suptitle('Continuous geometry constraints; IMP and fixed wall unchanged',fontsize=15)
    fig.savefig(folder/'geometry_constraints.png',dpi=180);plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    peak=[rows[0]['baseline_peak_MW_m2']]+[r['peak_MW_m2'] for r in rows]
    bars=axes[0].bar(['Baseline','20 mm','50 mm','100 mm'],peak,color=['#777777']+list(COLORS.values()))
    axes[0].bar_label(bars,fmt='%.3f',padding=3);axes[0].set(ylim=(0,max(peak)*1.18),ylabel='Peak heat flux [MW/m²]',title='Refined global wall peak')
    surface=next(s for s in study['metadata']['surfaces'] if s['name']=='floor')
    j=next(i for i,s in enumerate(study['metadata']['surfaces']) if s['name']=='floor');start=surface['segment_index']+j*(641-2)
    for mm,a in [(None,base)]+list(data.items()):
        wall,q=cells(a);mid=(wall[:-1]+wall[1:])/2;mask=(a['cell_parent']>=start)&(a['cell_parent']<start+640)
        axes[1].plot(mid[mask,0],q[mask],'--' if mm is None else '-',color='#555555' if mm is None else COLORS[mm],
            label='Baseline' if mm is None else f'{mm} mm')
    axes[1].set(xlabel='Floor R [m]',ylabel='Cell-average q [MW/m²]',title='Heat-flux distribution along floor');axes[1].legend(fontsize=9)
    for mm in COLORS:
        history=[json.loads(x) for x in (folder/f'R{mm:03d}'/'evaluations.jsonl').read_text().splitlines()]
        best=np.inf;x=[];y=[]
        for row in history:
            if row['feasible']:best=min(best,row['objective'])
            if np.isfinite(best):x.append(row['case']);y.append(best)
        axes[2].plot(x,y,color=COLORS[mm],label=f'{mm} mm')
    axes[2].axhline(1,color='#555555',ls='--');axes[2].set(xlabel='Unique candidate',ylabel='Best peak / baseline peak',title='Search history (320 tubes)');axes[2].legend(fontsize=9)
    for ax in axes:ax.grid(axis='y',alpha=.15)
    fig.savefig(folder/'peak_and_search.png',dpi=180);plt.close(fig)
    lines=['# OMP 접선·곡률·개구부 제약 재최적화','',
        '기준 D3D 형상은 원본 그대로 보존했습니다. OMP 왼쪽 바닥 접합부는 위치와 접선 방향을 고정하고, 가변 OMP 곡선에 절대 곡률 제한을 적용했습니다.',
        '오른쪽 고정 lip 끝점에서 바닥까지의 최소 거리는 기준 33.354336 mm 이상입니다. IMP 영역과 고정 벽은 변형하지 않았습니다.','',
        '| 최소 반경 제약 | peak (MW/m²) | 기준 대비 감소율 | 실제 바닥 최소 반경 | 최소 간격 | 격자 검증 | 개선 확인 |',
        '|---:|---:|---:|---:|---:|---|---|']
    for r in rows:lines.append(f"| {r['radius_constraint_mm']} mm | {r['peak_MW_m2']:.6f} | {r['reduction_percent']:.3f}% | {r['actual_floor_minimum_radius_mm']:.3f} mm | {r['minimum_gap_mm']:.3f} mm | {r['resolution_screen_passed']} | {r['improvement_confirmed']} |")
    lines += ['',f"원본 기준 peak: **{rows[0]['baseline_peak_MW_m2']:.6f} MW/m²**. 음의 감소율은 기준보다 peak가 증가했다는 뜻입니다.",
        '표 수치는 cubic 1280 tubes, geometry 641 nodes/surface, receiving wall 0.25 mm입니다. 검증은 320/640/1280 tubes, 별도 wall refinement, quintic 640을 포함합니다.',
        '원래 기준 형상은 새 접선 조건을 만족하는 후보가 아닙니다. 비교 기준과 후보 집합을 분리했으며, 조건을 어긴 기준 형상을 최적해로 반환하지 않습니다.','',
        '![열유속 비교](heat_flux_comparison.png)','', '![형상 제약 확대](geometry_constraints.png)','', '![peak와 탐색](peak_and_search.png)','',
        '## 고정 입력과 탐색 범위','',
        'OMP 열원 1 MW, n=1e19 m⁻³, χ⊥=0.3 m²/s, λq=3 mm, γ=7, κ₀=2000. 이전과 같은 단방향 열전도·sheath 폭 모델과 전력 보존형 벽 경계 분포를 사용했습니다.',
        'OMP의 네 구간 각각 4개 자유 계수를 탐색했습니다. 원래의 길이 1.2배 상한, 구조물·LCFS·자기 교차 금지, shelf 입사각 조건, ±100 mm 계수 범위를 유지했습니다.',
        '위·아래 양 방향 초기 형상과 확률 탐색, 좌표 미세 탐색을 사용했습니다. 실제 최고 성능은 기록된 탐색 예산에서 찾은 결과이며 전역 최적임을 증명하지 않습니다.',
        'GPU는 Gaussian 합산에 사용했고, 전력 CDF·형상 교차·열 상태·최대점 정밀 탐색은 기존 알고리즘을 유지했습니다. 최초 유효 후보와 최종 후보의 CPU/GPU peak 일치 여부를 기록했습니다.','',
        '## 해석 범위','',
        '접선 연결은 G1입니다. 접합부의 접선 방향은 연속이지만 곡률까지 연속인 G2는 강제하지 않았습니다. 요청한 왼쪽 접합부 외의 원래 고정 구조물 코너는 그대로입니다.',
        '간격은 고정 lip 끝점과 바닥 사이의 유클리드 최소 거리입니다. 개구부 전체의 유동·중성입자 배기 성능을 대체하는 모델은 아닙니다.',
        '최대 곡률과 최소 간격은 각 cubic span의 정상점과 끝점에서 검사했습니다. 교차 및 수송 도착 계산은 여전히 곡선을 근사한 polyline을 사용하므로 geometry node refinement를 함께 검증했습니다.',
        '수치 검증은 동일한 축약 수송 모델 안의 비교 신뢰도를 확인합니다. 실험값 일치나 완전한 플라즈마 수송 모델 검증을 의미하지 않습니다.','',
        '## 재현','',
        '`run_tangent_curvature_study.py --output <새 D 드라이브 경로>`로 탐색하고, `qualify_tangent_candidates.py --output <경로>`로 상위 후보를 640 tubes에서 다시 비교합니다. 이후 `verification/validate_tangent_curvature.py --output <경로>` 및 `visualize_tangent_curvature.py --output <경로>`를 실행합니다.',
        '640 tubes 후보 선별에는 각 탐색의 상위 5개 후보와 더 엄격한 반경 조건의 선택 후보를 포함합니다. 이 선별 결과와 원래 탐색 승자는 별도 기록합니다.',
        '각 R020/R050/R100 폴더의 manifest, evaluations.jsonl, search.json, validation.json과 refined CSV/NPZ에서 계수·형상·열유속·전력 장부·검증 결과를 확인할 수 있습니다.']
    (folder/'REPORT.ko.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
