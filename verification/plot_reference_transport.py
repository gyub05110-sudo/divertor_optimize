"""Reference-point report; diagnostics and failed screens remain explicit."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import LineString
from d3d_optimization.broadening_study import load_json


def data(folder,mode,candidate):
    for suffix in ['_refined.npz','.npz']:
        p=folder/mode/(candidate+suffix)
        if p.exists():return np.load(p)
    return None


def figure(folder,manifest):
    base=data(folder,'direct','baseline')
    if base is None:return
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    wall_ref=LineString(base['wall'])
    selections=[('direct','baseline','Direct baseline','#79828c','--'),
                ('direct','direct_best','Direct selected','#244e78','-'),
                ('thermal','baseline','Thermal baseline','#d6a461','--'),
                ('thermal','thermal_best','Thermal selected','#bd4e23','-'),
                ('frozen_transport','frozen_transport_best','Frozen-width selected','#6d5792',':')]
    for mode,candidate,label,color,style in selections:
        arr=data(folder,mode,candidate)
        if arr is None:continue
        wall=arr['wall'];p=arr['points'];q=arr['q_MW_m2'][:,0]
        axes[0,0].plot(wall[:,0],wall[:,1],ls=style,color=color,label=label)
        keep=np.isfinite(p).all(axis=1)&np.isfinite(q)&(q>max(q)*1e-5)&(p[:,1]<-1)
        s=shapely.line_locate_point(wall_ref,shapely.points(p[keep]));order=np.argsort(s)
        axes[0,1].plot(s[order],q[keep][order],ls=style,color=color,label=label)
        if 'sigma_m' in arr:
            x=arr['source_x_m']*1000;near=x<15
            axes[1,0].plot(x[near],arr['sigma_m'][near]*1000,ls=style,color=color,label=label)
            axes[1,1].plot(x[near],arr['thermal_Tt_eV'][near],ls=style,color=color,label=label)
    axes[0,0].set(xlim=(.98,1.85),ylim=(-1.45,-1.08),xlabel='R (m)',ylabel='Z (m)',title='Plate geometry')
    axes[0,0].set_aspect('equal',adjustable='box')
    axes[0,1].set(xlabel='Reference wall arc (m)',ylabel='Wall heat flux (MW/m²)',title='Absolute heat load')
    axes[1,0].set(xlabel='Upstream radial distance (mm)',ylabel='Deposited sigma in flux coordinate (mm)',title='Actual kernel widths')
    axes[1,1].set(xlabel='Upstream radial distance (mm)',ylabel='Pre-redistribution Tt (eV)',title='Tube closure temperature; no deposition feedback')
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0,1].legend(fontsize=8)
    fig.suptitle(f"{manifest['source']}: hypothetical reference; validity and regime diagnostics in report")
    fig.tight_layout();fig.savefig(folder/'comparison.png',dpi=160);plt.close(fig)


def write_report(output):
    output=Path(output)
    lines=['# 고정 기준 운전점의 형상 비교','',
        '입력 전력·상류 밀도·수송계수·sheath 규칙을 모든 형상에 동일하게 적용했다. 전도 two-point 온도 계산 후 횡방향 Gaussian 재분배를 하는 **one-way 근사**이며, 최종 열분포가 온도에 되먹임되는 2차원 수송 계산은 아니다. 실험 보정과 전역 최적해 증명은 하지 않았다.','',
        'direct는 기존 직접 매핑, frozen_transport는 기준 형상의 경로별 확산 폭을 고정한 대조군, thermal은 각 형상에서 온도와 확산 폭을 다시 계산하는 모델이다. 이번 목적함수는 **상류 폭 3 mm 한 조건의 피크**이며 과거 다섯 폭의 minimax 결과와 직접 비교하지 않는다.','',
        'IMP와 OMP는 각각 지정 전력을 받는 독립 시나리오다. 두 결과를 합쳐 실제 장치 총전력이라고 해석하지 않는다. 온도는 재분배 전 경로 진단값이며 최종 벽 위치의 플라즈마 온도가 아니다.','']
    for source in ['IMP','OMP']:
        folder=output/source
        if not (folder/'manifest.json').exists():continue
        m=load_json(folder/'manifest.json');ref=m['reference']
        v=load_json(folder/'validation.json') if (folder/'validation.json').exists() else {}
        status='검증 완료' if v.get('complete') else ('검증 진행 중' if v else '검증 전')
        lines += [f'## {source}', '',f"기준: P={ref['power_W']/1e6:g} MW, nu={ref['density_m3']:.3g} m⁻³, χ⊥={ref['chi_perp_m2_s']:g} m²/s, λq={ref['lambda_q_m']*1000:g} mm. {status}. 탐색 {m['metadata']['tubes']} 경로; smoke={m['smoke_test']}.",'',
            '| 모델 | 탐색 기준 적합 | 기준 피크 MW/m² | 탐색 미해결 전력 | 선택 후보/기준 | 수치 판정 |',
            '|---|---|---:|---:|---:|---|']
        for mode,model in m['models'].items():
            base=model['baseline'];candidate=mode+'_best'
            raw=m.get('cross_evaluation',{}).get(mode,{}).get(candidate,{})
            verdict=v.get('comparisons',{}).get(mode,{}).get(candidate,{})
            ratio=verdict.get('finest_peak_ratio',raw.get('worst_peak_ratio'))
            peak=base.get('peak_q_MW_m2',[None])[0]
            missing=base.get('missing_power')
            fmt=lambda x:'—' if x is None else f'{x:.6g}'
            decision='수치 개선 확인' if verdict.get('improvement_confirmed') else verdict.get('verification_status','검증 전')
            if not base['feasible']:decision='기준 부적합: 탐색 제외'
            lines.append(f"| {mode} | {base['feasible']} | {fmt(peak)} | {fmt(None if missing is None else 100*missing)}% | {fmt(ratio)} | {decision} |")
        lines+=['','비율은 같은 모델·같은 운전점의 기준 형상 대비다. 검증 비율이 없으면 탐색 비율을 표시하며, 실패한 행의 비율은 개선 증거가 아니다. 전력 분배 영역 밖/차폐 전력이 0.1%를 초과하면 부적합이며 재정규화하지 않는다.','']
        lines+=['### 동일 후보의 모델 간 교차 평가','',
                '| 후보 | direct | frozen_transport | thermal |',
                '|---|---:|---:|---:|']
        refined=v.get('complete') and len(v.get('levels',[]))>=3
        cross=v['levels'][2]['models'] if refined else m.get('cross_evaluation',{})
        for candidate in m['candidates']:
            values=[]
            for mode in ['direct','frozen_transport','thermal']:
                row=cross.get(mode,{}).get(candidate,{})
                ratio=row.get('worst_peak_ratio')
                values.append('—' if ratio is None else f"{ratio:.6f} ({'적합' if row.get('feasible') else '부적합'})")
            lines.append('| '+candidate+' | '+' | '.join(values)+' |')
        lines+=['',f"{'가장 미세한 cubic 검증' if refined else '탐색'} 격자의 후보/동일 모델 기준 피크 비율이다. 적합 표시는 해당 해상도에서의 제약 판정이며, 전 격자 수치 개선 확인은 위 판정을 따른다.",'']
        for mode in m['models']:
            rows=[lev['models'].get(mode,{}).get('baseline',{}) for lev in v.get('levels',[])]
            misses=[r['missing_power'] for r in rows if 'missing_power' in r]
            if misses:lines += [f"- {mode}: 검증 기준 최대 미해결 전력 {100*max(misses):.6g}%; 모든 기준 audit 통과={all(r.get('audit_passed',False) for r in rows)}."]
        lines+=['','### 수송 근사의 적용 상태','']
        thermal=m['models'].get('thermal',{}).get('baseline',{}).get('thermal',{})
        if thermal:
            lines += [f"기준 형상의 전력 가중 Tu={thermal['Tu_eV']['power_weighted_mean']:.3f} eV, 재분배 전 Tt={thermal['Tt_eV']['power_weighted_mean']:.3f} eV, sigma={thermal['sigma_m']['power_weighted_mean']*1000:.4f} mm.",'',
                f"Tt<5 eV 경로의 전력 비중={100*thermal['cold_target_power_fraction']:.3f}%; 전도 온도강하 지표<0.1인 경로 비중={100*thermal['weak_conduction_power_fraction']:.3f}%; sigma>0.3λq 비중={100*thermal['sigma_above_0p3_lambda_power_fraction']:.3f}%.",'',
                f"벌크 free-streaming 열유속 대비 비율>0.1인 경로의 전력 비중={100*thermal['free_streaming_ratio_above_0p1_power_fraction']:.3f}%. 이 진단은 전자 유동이나 운동론 해석을 대체하지 않는다.",'',
                '낮은 Tt는 이 모델의 적용 한계 경고이며 detachment 판정이 아니다. 전도 온도강하가 작다고 유동을 해석한 것이 아니며, 전도 지배의 타당성은 별도 검증이 필요하다. 극저전력 꼬리의 비현실적 상태도 자르거나 온도 하한으로 숨기지 않았다.','']
        if v.get('sensitivity'):
            lines+=['### 고정 후보의 운전점·국소 분배 가정 민감도','',
                    '| 조건 | 후보 | 기준 적합 | 후보 적합 | 후보/기준 피크 |','|---|---|---|---|---:|']
            for label,rec in v['sensitivity'].items():
                for candidate,row in rec['candidates'].items():
                    ratio=row.get('worst_peak_ratio')
                    lines += [f"| {label} | {candidate} | {rec['baseline']['feasible']} | {row['feasible']} | {'—' if ratio is None else f'{ratio:.6f}'} |"]
            lines+=['','이 표는 동일 후보를 한 해상도에서 재평가한 민감도이며, 재최적화나 전 해상도 검증 결과가 아니다. 부적합 기준 또는 후보의 비율은 진단용이다.','']
        figure(folder,m)
        if (folder/'comparison.png').exists():lines += [f'![{source}]({source}/comparison.png)','']
        lines += [f'원자료: [{source} manifest]({source}/manifest.json), [{source} validation]({source}/validation.json).','']
    lines+=['## 재현과 한계','',
        '설정은 각 manifest, 계산 전 코드·입력 해시는 provenance.json에 보존한다. 기존 direct 및 길이 의존 broadening 코드는 변경하지 않았다. 수치 검증은 first hit, 독립 자기장/열적분, 전력 장부/면적 적분, 형상 제약, 격자·보간법 민감도를 확인한다.','',
        '정적 bulk pressure와 sonic sheath의 factor-two 관계를 쓰며 입자·운동량·중성입자·복사·분리(detachment)를 풀지 않는다. 최종 재분배 열유속은 재분배 전 Tt의 sheath 식을 다시 만족시키도록 계산되지 않는다. 자기장 고정 plate 비교이며 Super-X 자기배위 설계가 아니다.','',
        'Spreading은 [Nille–Sieglin–Eich 전도 근사](https://arxiv.org/html/1610.04148)의 적분을 전 경로/가변 자기장에 확장한 가정이다. χ⊥는 전자 열용량 1.5 e n 기준이다.','']
    output.mkdir(parents=True,exist_ok=True)
    (output/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    write_report(parser.parse_args().output)
