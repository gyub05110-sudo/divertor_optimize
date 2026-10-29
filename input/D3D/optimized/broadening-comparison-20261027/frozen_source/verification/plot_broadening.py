"""Static, inspectable comparison figures and a scoped scientific report."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import LineString
from d3d_optimization.broadening_study import load_json


def arrays(folder, model, candidate):
    refined = folder/model/(candidate+'_refined.npz')
    raw = folder/model/(candidate+'.npz')
    return np.load(refined if refined.exists() else raw) if refined.exists() or raw.exists() else None


def source_figure(folder, manifest):
    available = [(name, model) for name,model in manifest['models'].items()
                 if model['config']['sigma_ref_m'] > 0 and name+'_best' in manifest['candidates']]
    if not available:
        return None
    model = available[0][0]
    direct_name = 'direct_best' if 'direct_best' in manifest['candidates'] else 'historical_direct'
    selections = [('direct','baseline','Direct: baseline','#7f8c8d','--'),
                  ('direct',direct_name,'Direct: selected','#183d68','-'),
                  (model,'baseline',f'{model}: baseline','#d9a760','--'),
                  (model,model+'_best',f'{model}: selected','#bd4e23','-')]
    fig,axes = plt.subplots(1,3,figsize=(15,4.7),gridspec_kw={'width_ratios':[1,1.3,1.1]})
    reference = arrays(folder,'direct','baseline')
    if reference is None:
        plt.close(fig)
        return None
    line = LineString(reference['wall'])
    for current,candidate,label,color,style in selections:
        data = arrays(folder,current,candidate)
        if data is None:
            continue
        wall = data['wall']
        axes[0].plot(wall[:,0],wall[:,1],color=color,ls=style,lw=1.6,label=label)
        points,heat = data['points'],data['q_over_P'][:,2]
        keep = np.isfinite(points).all(axis=1)&np.isfinite(heat)&(heat > max(heat)*1e-5)&(points[:,1]<-1.)
        x = shapely.line_locate_point(line,shapely.points(points[keep]))
        order = np.argsort(x)
        axes[1].plot(x[order],heat[keep][order],color=color,ls=style,lw=1.4,label=label)
        axes[2].plot([1,2,3,5,10],data['peak_q'],color=color,ls=style,marker='o',label=label)
    axes[0].set(xlabel='R (m)',ylabel='Z (m)',title='Plate geometry',xlim=(.98,1.85),ylim=(-1.45,-1.08))
    axes[0].set_aspect('equal',adjustable='box')
    axes[1].set(xlabel='Reference-wall arc coordinate (m)',ylabel='q / P (1/m²)',title='Heat profile, upstream λq = 3 mm')
    axes[2].set(xlabel='Upstream λq (mm)',ylabel='All-wall peak q / P (1/m²)',title='Peaks at five source widths')
    for ax in axes:
        ax.grid(alpha=.2)
    axes[1].legend(fontsize=7)
    fig.suptitle(f"{manifest['source']} — exploratory deposition surrogate; see numerical validity in report",fontsize=12)
    fig.tight_layout()
    filename = folder/'comparison.png'
    fig.savefig(filename,dpi=170)
    plt.close(fig)
    return filename


def write_report(output):
    output = Path(output)
    lines = ['# 연결 길이 의존 broadening 비교', '',
        '고정 자기장과 기존 24개 벽 형상 자유도에서 수행한 탐색이다. 새 모델은 열침적 근사이며 평행·횡방향 수송 방정식, 입자 수송, 복사 또는 detachment를 푼 결과가 아니다.', '',
        '기존 지수 열원과 단위 입력 전력을 유지한다. 실제 3차원 source-to-first-hit 연결 길이로 Gaussian 표준편차를 정하고, 같은 단조 자기선속 벽 가지에 분배한다. 전력은 전벽 침적·차폐·가지 밖·미충돌로 나누며 남은 전력을 재정규화하지 않는다.', '',
        '`sigma_ref`는 upstream-equivalent 좌표의 표준편차이며 Eich S나 측정 확산계수가 아니다. 기본 L_ref=10 m, 지수=1, anchor=5 mm poloidal이다. 아래 수치는 물리적 최적해나 Super-X의 성능 예측을 뜻하지 않는다.', '']
    for source in ['IMP','OMP']:
        folder = output/source
        if not (folder/'manifest.json').exists():
            continue
        manifest = load_json(folder/'manifest.json')
        validation = load_json(folder/'validation.json') if (folder/'validation.json').exists() else None
        lines += [f'## {source}', '',
            f"모드: {manifest['mode']}; smoke: {manifest['smoke_test']}; 탐색 해상도: {manifest['metadata']['tubes']} 경로. 두 열원은 각각 단위 전력이며 합산하지 않는다.", '',
            '| 모델 | 기준 형상 적합 | 미해결 전력 최대 | 선택 후보 / 기준 최대 비율 | 수치 개선 확인 |',
            '|---|---:|---:|---:|---|']
        for model,record in manifest['models'].items():
            baseline = record['baseline']
            candidate = model+'_best'
            row = manifest.get('cross_evaluation',{}).get(model,{}).get(candidate,{})
            verdict = validation.get('comparisons',{}).get(model,{}).get(candidate,{}) if validation else {}
            ratio = verdict.get('finest_worst_peak_ratio',row.get('worst_peak_ratio'))
            missing = baseline.get('missing_power')
            ratio_text = '—' if ratio is None else f'{ratio:.6f}'
            missing_text = '—' if missing is None else f'{100*missing:.6g}%'
            status = ('확인' if verdict.get('improvement_confirmed') else verdict.get('verification_status','검증 전'))
            if not baseline['baseline_admissible']:
                status = '기준 부적합: 탐색 제외'
            lines.append(f"| {model} | {baseline['baseline_admissible']} | {missing_text} | {ratio_text} | {status} |")
        lines += ['', '각 비율은 **같은 모델의 기준 형상**과 비교한 다섯 λq의 최대값이다. 미해결 전력이 0.1%를 초과하면 해당 모델/형상은 부적합이며, 기준 형상이 부적합하면 개선율을 인정하지 않는다.', '']
        for model,record in manifest['models'].items():
            selected = model+'_best'
            if selected not in manifest['candidates'] or record['config']['sigma_ref_m']==0:
                continue
            cross = manifest['cross_evaluation'][model]
            choice = 'direct_best' if 'direct_best' in cross else 'historical_direct'
            raw_selected = cross[selected]
            raw_direct = cross.get(choice,{})
            lines += [f'### {model}: 형상 교차 평가', '',
                f"탐색 격자에서 새 후보의 기준 대비 최대 형상 변위: {1000*raw_selected.get('max_wall_displacement_m',0):.4f} mm.", '',
                f"직접 매핑 모델의 선택 후보와의 최대 형상 차이: {1000*raw_selected.get('max_wall_difference_from_direct_selected_m',0):.4f} mm.", '',
                '| 같은 broadening 모델에서 평가 | 최대 q/P 비율 | 후보 적합 |',
                '|---|---:|---|']
            for name,row in [(choice,raw_direct),(selected,raw_selected)]:
                ratio=row.get('worst_peak_ratio')
                lines.append(f"| {name} | {'—' if ratio is None else format(ratio,'.6f')} | {row.get('feasible',False)} |")
            if validation:
                v=validation['comparisons'][model].get(selected,{})
                lines += ['',f"세부 검증: `{v.get('verification_status')}`; 관측 비율 변화={v.get('max_ratio_change')}; 절대 피크 상대 변화={v.get('max_absolute_peak_relative_change')}; anchor 민감={v.get('anchor_sensitive')}. 수치 오차와 모델 가정 민감도를 구분한다."]
            lines += ['']
        figure=source_figure(folder,manifest)
        if figure:
            lines += [f'![{source} comparison]({source}/comparison.png)', '']
        lines += [f'상세 자료: [{source} manifest]({source}/manifest.json), 원시 탐색 로그는 각 모델 폴더의 `evaluations.jsonl`.', '']
        if validation:
            lines += [f'전체 격자·보간법·차폐거리·길이 지수 비교: [{source} validation]({source}/validation.json).', '']
    lines += ['## 재현 및 해석', '',
        '실제 실행 인수는 각 manifest의 settings, 코드와 입력 SHA256은 provenance.json에 저장한다. 검증기는 변경된 코드/입력으로 과거 탐색을 검증했다고 주장하지 않도록 해시 불일치를 거부한다.', '',
        '연결 길이를 늘렸을 때 σ가 늘어나도록 정의한 현상론적 모델이다. 길이에 대한 선호가 나타나도 그것만으로 평행 수송의 물리적 예측력이 검증되지는 않는다. 같은 자기장에서의 plate 비교이며 자기배위 최적화가 아니다.', '',
        '전력 장부 보존과 독립 면적 적분, 격자 수렴은 서로 다른 검사다. 수치 검사 통과는 실험 검증이나 전역 최적해 증명이 아니다. 강한 퍼짐에서 가지 밖 전력이 커지면 더 넓은 물리적 수송 모델이 필요하다.', '']
    output.mkdir(parents=True,exist_ok=True)
    (output/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    write_report(p.parse_args().output)
