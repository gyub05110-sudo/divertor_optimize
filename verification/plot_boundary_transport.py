"""Static scientific plots and Korean result report from saved evidence."""
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from d3d_optimization.broadening_study import load_json


def report(folder,benchmark=None):
    folder=Path(folder);m=load_json(folder/'manifest.json');v=load_json(folder/'validation.json')
    base=np.load(folder/'baseline_refined.npz');selected=np.load(folder/'selected_refined.npz')
    rows=v['levels'][2]['rows'];colors=['#2876b5','#db682a']
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for name,a,color in [('baseline',base,colors[0]),('selected',selected,colors[1])]:
        wall=a['wall'];axes[0,0].plot(wall[:,0],wall[:,1],color=color,label=name,lw=1.6)
        receive=a['receiving_wall'];ds=np.linalg.norm(np.diff(receive,axis=0),axis=1)
        arc=np.r_[0.,np.cumsum(ds)];mid=(arc[:-1]+arc[1:])/2
        q=a['cell_power']/np.diff(a['area_edges'])*m['reference']['power_W']/1e6
        axes[0,1].plot(mid,q,color=color,label=name)
        mask=q>max(rows['baseline']['peak_MW_m2'],rows['selected']['peak_MW_m2'])*.01
        if np.any(mask):
            bounds=(mid[mask].min()-.05,mid[mask].max()+.05)
            if name=='baseline':qb=bounds
            else:axes[0,1].set_xlim(min(qb[0],bounds[0]),max(qb[1],bounds[1]))
        axes[1,1].plot(a['source_x_m']*1000,a['sigma_physical_m']*1000,color=color,label=name)
    axes[0,0].set(xlim=(.95,1.65),ylim=(-1.45,-.95),xlabel='R [m]',ylabel='Z [m]',title='Fixed baseline and selected divertor geometry')
    axes[0,0].set_aspect('equal');axes[0,0].legend()
    axes[0,1].set(xlabel='Connected wall arclength [m]',ylabel='Cell-average heat flux [MW/m²]',title='Arriving heat profile (1280 donors)');axes[0,1].legend()
    names=list(rows['baseline']['region_power']);x=np.arange(len(names))
    for j,n in enumerate(['baseline','selected']):
        axes[1,0].bar(x+(j-.5)*.36,[100*rows[n]['region_power'][k] for k in names],.36,label=n,color=colors[j])
    axes[1,0].set_xticks(x,names,rotation=22,ha='right');axes[1,0].set(yscale='log',ylim=(1e-7,120),ylabel='Input power [%]',title='Power received by each material region');axes[1,0].legend()
    axes[1,1].set(xlim=(0,20),yscale='log',xlabel='OMP source offset [mm]',ylabel='Effective physical footprint σ [mm]',title='Geometry-dependent width, identical operating inputs')
    visible=np.concatenate([a['sigma_physical_m'][a['source_x_m']<=.02]*1000 for a in [base,selected]])
    visible=visible[np.isfinite(visible)&(visible>0)]
    axes[1,1].set_ylim(visible.min()*.8,visible.max()*1.25)
    for ax in axes.ravel():ax.grid(alpha=.2)
    fig.suptitle(f"OMP | {rows['baseline']['peak_MW_m2']:.4f} → {rows['selected']['peak_MW_m2']:.4f} MW/m² | area-boundary footprint model")
    fig.savefig(folder/'comparison.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,6),layout='constrained',sharex=True,sharey=True)
    vmax=rows['baseline']['peak_MW_m2']
    for ax,a,label in zip(axes,[base,selected],['Fixed baseline','Selected geometry']):
        w=a['receiving_wall'];q=a['cell_power']/np.diff(a['area_edges'])*m['reference']['power_W']/1e6
        ax.plot(a['wall'][:,0],a['wall'][:,1],color='.8',lw=.8)
        lc=LineCollection(np.stack([w[:-1],w[1:]],axis=1),cmap='inferno',norm=plt.Normalize(0,vmax),linewidth=4)
        lc.set_array(q);ax.add_collection(lc);ax.set(xlim=(.95,1.7),ylim=(-1.45,-.9),title=label,xlabel='R [m]',ylabel='Z [m]');ax.set_aspect('equal')
    fig.colorbar(lc,ax=axes,label='Cell-average heat flux [MW/m²]');fig.savefig(folder/'arrival_map.png',dpi=180);plt.close(fig)
    lines=['# 고정 D3D 기준 형상에서 OMP 최적화','',
      '목표는 주어진 형상과 가상 운전 입력을 사용하는 빠른 형상 비교·최적화입니다. 실제 측정값과의 일치는 목표와 판정 기준에 포함하지 않았습니다.','',
      '## 결과','',f"- 기준 최대 열유속: **{rows['baseline']['peak_MW_m2']:.6f} MW/m²**",
      f"- 선택 형상 최대 열유속: **{rows['selected']['peak_MW_m2']:.6f} MW/m²**",
      f"- 감소율: **{v['verdict'].get('finest_reduction_percent',float('nan')):.3f}%**",
      f"- 수치 검증에서 개선 확인: **{v['verdict']['improvement_confirmed']}**",
      f"- 탐색: {m['search']['evaluated']}개 고유 후보, {m['search']['feasible_evaluations']}개 적합, {m['search']['elapsed_seconds']:.1f}초",
      f"- 전체 탐색 실행(준비·저장 포함): {m['elapsed_seconds']:.1f}초",
      f"- 기준 형상 준비: {m['setup_seconds']:.1f}초, 기준 형상 1회 평가: {m['baseline']['elapsed_seconds']:.2f}초",
      f"- 최대 형상 변위: {rows['selected']['max_wall_displacement_m']*1000:.3f} mm",
      f"- 전력가중 평균 연결 길이: {rows['baseline']['thermal']['length_m']['power_weighted_mean']:.3f} → {rows['selected']['thermal']['length_m']['power_weighted_mean']:.3f} m",
      '- 전역 최적해를 증명한 결과가 아니라 주어진 탐색 예산에서 찾은 적합 후보입니다.','',
      '## 입력과 해석 범위','',f"`{m['reference']}`",'',
      f"P는 선택한 {m['source']} 소스의 계산 구간(0–{m['metadata']['support_mm']:.3f} mm)에 부여한 전력입니다. 전체 장치 전력이나 IMP/OMP 동시 전력 분할을 뜻하지 않습니다.",
      '기준 벽과 자기 평형은 보존합니다. 각 후보에서 첫 도착점·연결 길이·열적 폭을 다시 계산합니다. Te는 지정한 P, n, sheath/전도 계수에 따른 two-point closure의 출력이며 독립적으로 강제한 입력은 아닙니다.','',
      '## 무엇을 바꿨는가','',
      '기존 OMP320 미배분 0.927699%는 자속 구간 밖 꼬리 0.796088%와 오목한 벽에서의 직선 차폐 제외 0.131611%였습니다. 기준 형상이 잘못된 것이 아니라 기존 받는 영역의 근사에서 제외된 전력입니다.',
      '새 모델은 첫 도착점을 중심으로 연결된 닫힌 벽의 누적 회전면적 A=∫2πR ds에 주기 Gaussian을 정의합니다. 각 셀의 전력은 CDF 차이로 계산하여 도착한 전력을 재정규화 없이 보존합니다. 첫 벽에 도착하지 않은 전력은 미도착으로 남깁니다.',
      '이것은 명시적인 **벽 footprint 대리모델**입니다. 기존 투영 Gaussian의 단순 범위 확장이 아닙니다. σ_flux를 국소 |∇ψ|로 변환한 물리 폭을 유효 벽방향 폭으로 사용하며, 1/sin 입사각에 의한 폭 투영은 적용하지 않습니다. 따라서 실제 플라스마 속 횡방향 궤적을 복원했다고 해석하면 안 됩니다.','',
      '## 수치 검증','',
      '| 보간 / donor 수 | 기준 qmax | 선택 qmax | 선택/기준 | 독립 audit |','|---|---:|---:|---:|---|']
    for l in v['levels']:
        b,s=l['rows']['baseline'],l['rows']['selected']
        lines.append(f"| {l['method']} / {l['tubes']} | {b.get('peak_MW_m2',float('nan')):.6f} | {s.get('peak_MW_m2',float('nan')):.6f} | {s.get('peak_ratio',float('nan')):.7f} | {b['audit_passed'] and s['audit_passed']} |")
    lines.extend(['','전력 장부와 면적 적분, 독립 첫 교차점, 독립 자기장·열적 적분, 고정 벽 좌표, 원래의 형상/입사각 제약, 연속 spline 길이를 확인했습니다. 320/640/1280 donor 및 quintic 보간, 받는 벽만 세분화한 별도 검증을 포함합니다.',
      '',f"판정 상세: `{v['verdict']}`",'',
      '## 남는 한계','',
      '- 고정 자기장, 단방향 전도/sheath closure입니다. 복사·재결합·중성입자·운동량·시간 변화 및 횡수송이 온도에 주는 피드백은 풀지 않습니다.',
      '- 폭이 매우 넓으면 연결된 전체 벽으로 전력이 퍼지는 모델입니다. 폭/전체 면적 비율과 0.1 m 이상 이동한 전력을 함께 기록합니다.',
      '- 지원 대상은 하나의 단순한 닫힌 양의 R 벽입니다. 서로 분리된 벽 성분 사이의 전달은 구현하지 않았습니다.',
      '- 최대점 탐색은 donor 근처와 벽 격자의 샘플에 연속 최적화를 결합합니다. 수학적인 전역 최대점 증명은 아니므로 수치 검증을 유지해야 합니다.',
      '- 연결 길이가 길수록 무조건 좋은 목적함수를 넣지 않았습니다. 이번 선택 형상의 평균 연결 길이는 오히려 짧아졌습니다. 열유속 감소에는 기하·도착점·폭 변환이 함께 작용하므로 감소율을 parallel transport 단독 효과나 Super-X의 우월성으로 해석할 수 없습니다.',
      '- 양의 유한 입력을 받지만 float 분해능 이하의 폭, overflow, 0 자기장 기울기 등은 명시적으로 실패합니다. 모든 수학적 입력에서 적합 형상이나 개선이 존재한다고 보장하지 않습니다.',
      '- 평균 벽 이동 거리는 셀 중심 근사이며 최대 셀 길이/2의 오차 상한을 기록합니다. 0.1 m 이상 이동한 전력은 CDF 구간으로 계산합니다.','',
      '![geometry and profiles](comparison.png)','', '![arrival heat map](arrival_map.png)','',
      '실행 방법은 저장소의 `BOUNDARY_OPTIMIZATION.md`, 재현 입력·계수는 이 폴더의 `manifest.json`, 검증은 `validation.json`, 상세 배열은 `*_refined.npz`에 있습니다.'])
    if benchmark:
        b=load_json(Path(benchmark)/'benchmark.json');lines.extend(['','## 기존 모델과의 고정 기준 형상 비교','',
            '| 소스 | 기존 미배분 [%] | 새 미배분 [%] | 기존 qmax [MW/m²] | 새 qmax [MW/m²] |','|---|---:|---:|---:|---:|'])
        for source,r in b['models'].items():
            lines.append(f"| {source} | {100*r['legacy']['missing_power']:.6f} | {100*r['new']['ledger']['unresolved']:.6f} | {r['legacy']['peak_q_MW_m2'][0]:.6f} | {r['new']['peak_MW_m2']:.6f} |")
        lines.extend(['','동일 형상·입력의 비교이며 새 모델은 전체 footprint 정의가 바뀌었습니다. 기존 미배분만 특정 표면에 더한 결과가 아닙니다.','', '## 입력 변화 점검','',
            '각 입력을 따로 바꾼 고정 형상 시험입니다. 아래 후보는 기존 historical 형상이며 이 연구의 선택 형상을 재최적화한 결과는 아닙니다.',
            '', '| 입력 변경 | 기준 적합 | 기존 후보 적합 | 후보/기준 qmax |','|---|---|---|---:|'])
        for k,r in b['input_sweep'].items():lines.append(f"| {k} | {r['baseline']['feasible']} | {r['historical_shape']['feasible']} | {r['historical_shape'].get('peak_ratio',float('nan')):.5f} |")
    path=folder/'REPORT.ko.md';path.write_text('\n'.join(lines)+'\n',encoding='utf-8');return path


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=Path,required=True);p.add_argument('--benchmark',type=Path)
    a=p.parse_args();print(report(a.folder,a.benchmark))
