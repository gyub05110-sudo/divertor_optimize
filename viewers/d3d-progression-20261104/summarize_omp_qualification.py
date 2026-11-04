"""Reports/figures from committed archived-candidate qualification only."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from qualify_omp_stages import read,sha,save
from export_qualified_finals import HERE,STUDY


def main():
    checked=read(HERE/'qualified-final-verification.json');index=read(HERE/'dist/qualified-finals.json');report=read(STUDY/'summary.json')
    if not checked['passed'] or checked['summary_sha256']!=sha(STUDY/'summary.json') or checked['overlay_sha256']!=sha(HERE/'dist/qualified-finals.json'):raise ValueError('verification changed')
    dest=STUDY/'report';dest.mkdir(exist_ok=True);data=read(HERE/'dist/data.json')
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');norm=Normalize(0,data['baseline']['refined']['OMP'])
    for mm,j in [(20,0),(50,1)]:
        if not any(e['final']['candidate']['radius_mm']==mm for e in index['entries']):
            for ax in axes[:,j]:
                ax.text(.5,.5,f'OMP Rmin {mm} mm\nNo archived candidate passed all existing criteria',ha='center',va='center',transform=ax.transAxes)
                ax.set_axis_off()
    for entry in index['entries']:
        j=[20,50].index(entry['final']['candidate']['radius_mm'])
        f=entry['final'];p=read(HERE/'dist'/entry['profile']['profile']);segments=np.array(p['segments']);wall=np.array(f['wall']);base=np.array(data['baseline']['wall'])
        ax=axes[0,j];ax.plot(base[:,0],base[:,1],color='.6',ls='--',lw=1,label='Fixed D3D baseline')
        ax.add_collection(LineCollection(segments[:,:4].reshape(-1,2,2),array=segments[:,4],cmap='inferno',norm=norm,linewidths=3))
        ax.set(xlim=(.97,1.73),ylim=(-1.44,-1.08),xlabel='R [m]',ylabel='Z [m]',title=f"OMP Rmin {f['candidate']['radius_mm']} mm | qmax {f['peak']:.6f} MW/m²\nFully qualified: {100*(1-f['ratio']):.3f}% reduction, archived stage {f['archived_stage']}")
        ax.set_aspect('equal');ax.grid(alpha=.15);ax.legend(loc='lower left',fontsize=8)
        profile=axes[1,j]
        for label,npz,color in [('Baseline',STUDY/'baseline/cubic_1280/baseline.npz','#777'),('Selected',Path(f['origin']).parent/'selected_refined.npz','#a12c56')]:
            with np.load(npz) as z:
                w=z['receiving_wall'];length=np.linalg.norm(np.diff(w,axis=0),axis=1);mask=(w[:-1,1]<-1.05)&(w[1:,1]<-1.05)
                centers=np.cumsum(length)-.5*length;centers-=centers[mask].min()
                q=z['cell_power']*data['reference']['power_W']/1e6/np.diff(z['area_edges'])
                profile.plot(centers[mask],q[mask],color=color,label=label,lw=1.6)
        profile.set(xlabel='Arc length along lower wall [m] (each shape)',ylabel='Cell-average q [MW/m²]');profile.grid(alpha=.15);profile.legend()
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='inferno'),ax=axes[0,:],label='Cell-average q [MW/m²]',shrink=.85)
    fig.suptitle('Additional qualification of archived 72-hour OMP milestones | Fixed 1 MW reference',fontsize=15)
    fig.savefig(dest/'omp-qualified-comparison.png',dpi=180);plt.close(fig)
    lines=['# OMP 보관 후보 추가 정밀검증 (2026-11-04)','',
        '새로운 형상 탐색 없이, 72시간 탐색에서 현재 형상 검사를 통과한 기록 갱신 후보 53개를 재검증했다. '
        '고정 D3D 기준 형상과 OMP 독립 1 MW 및 기존 수송 입력을 유지했다.','',
        '| Rmin | 전체 후보 | 선택 단계 | 정밀 peak (MW/m²) | 기준 대비 감소 | 정밀 순위 |',
        '|---:|---:|---:|---:|---:|---:|']
    for row in checked['cases']:
        lines.append(f"| {row['radius_mm']} mm | {row['ranked_candidates']} | {row['archived_stage']} | {row['peak_MW_m2']:.9f} | {row['reduction_percent']:.6f}% | {row['selected_fine_rank']} |")
    lines+=['',f"정밀 기준 peak: {data['baseline']['refined']['OMP']:.12f} MW/m².",'',
        '## 선택 방식과 검증','',
        '모든 53개 후보를 cubic 1280 경로에서 평가하여 peak 순서를 정했다. 각 반경에서 낮은 순서로 '
        '전체 검증을 수행했고 첫 통과 후보를 선택했다. 더 낮은 후보가 탈락했다면 그 근거는 summary.json의 attempts와 validation에 보존된다. '
        '아직 전체 검증을 하지 않은 나머지 후보는 선택 후보보다 1280 peak가 낮지 않다.','',
        '전체 검증은 cubic 320/640/1280, quintic 640, 별도 벽 격자 세분화, 독립 첫 도착점과 열적 적분, '
        '접선·곡률·간격·연속 벽 교차 검사, 반대 영역 보존, 전력 장부, CPU/GPU 열분포 및 셀 전력 비교를 포함한다. '
        '수치 검증 통과와 기준 대비 개선 확인을 모두 요구한다.','']
    lines+=['## 더 낮은 peak 후보의 탈락 이유','',
        '| Rmin | 전체 검증 시도 | 탈락 후보 | 선택된 정밀 순위 |','|---:|---:|---:|---:|']
    for mm in (20,50):
        attempts=report['attempts'][str(mm)];passed=[a for a in attempts if a['passed']]
        lines.append(f"| {mm} mm | {len(attempts)} | {sum(not a['passed'] for a in attempts)} | {passed[0]['fine_rank'] if passed else '통과 없음'} |")
    diagnostic=STUDY/'diagnostics/peak-location.json'
    if diagnostic.exists():
        diag=read(diagnostic);v=diag['verdict']
        lines+=['',
            f"20 mm 18단계는 peak 크기 변화 {100*v['max_absolute_peak_relative_change']:.4f}%와 누적 분포 변화 {v['max_wall_cdf_change']:.6f}는 허용 범위였지만, 최대점 위치 변화 {1000*v['max_peak_location_change_m']:.3f} mm가 기존 2 mm 기준을 넘었다. "
            '서로 떨어진 두 봉우리의 높이가 거의 같아 보간 방식에 따라 최대 봉우리가 바뀌는 현상이다. 전체 분포가 그 거리만큼 이동했다는 뜻이 아니다. '
            '형상·전력 보존 탈락과 구분해야 하며, 이번 작업에서는 최대점 위치 기준을 완화하지 않았다.','',
            '이 기준은 peak를 평탄화해 여러 봉우리 높이를 비슷하게 만드는 연구 목표와 충돌할 수 있다. '
            '후속 방법론 검토에서는 단일 argmax 위치와 near-maximum 영역·전력 누적 분포의 안정성을 구분할 필요가 있다. '
            '현재 최종 후보에는 기존 모든 통과 기준을 그대로 적용했다.','',
            '![20 mm 18단계의 근접 최대값 확대](../diagnostics/omp20-competing-peaks.png)','']
    for row in checked['cases']:
        v=row['verdict'];g=row['geometry']
        radius=min(r for r in g['actual_minimum_radius_m'].values() if r is not None)
        lines += [f"### Rmin {row['radius_mm']} mm",'',f"후보 ID: `{row['identity']}`",'',
            f"- 전체 검증 시도: {row['full_qualification_attempts']}개 / 1280 평가 {row['ranked_candidates']}개",
            f"- 실제 최소 곡률반경: {1000*radius:.6f} mm; 접선 오차: {g['tangent_angle_error_deg']:.3g}°",
            f"- 개구부 간격: {1000*g['gap_m']:.6f} mm (하한 {1000*g['gap_minimum_m']:.6f} mm)",
            f"- 반대 영역 최대 좌표 변화: {g['opposite_coordinate_error_m']:.3g} m",
            f"- P/qmax로 정의한 유효 면적: 기준 {1/row['baseline_MW_m2']:.6f} → 후보 {1/row['peak_MW_m2']:.6f} m² (단일 1 MW 입력)",
            f"- 원시 셀 전력 합: {row['raw_power_fraction']:.16g}; 화면 표시 셀 전력 합: {row['display_power_fraction']:.16g}",
            f"- 해상도 간 peak 최대 상대 변화: {100*v['max_absolute_peak_relative_change']:.6f}% (기준 0.5%)",
            f"- peak 비율 최대 변화: {v['max_peak_ratio_change']:.8g} (기준 0.002)",
            f"- 영역별 전력 L1 변화: {v['max_region_power_l1_change']:.8g} (기준 0.005)",
            f"- 누적 분포 변화: {v['max_wall_cdf_change']:.8g} (기준 0.002)",
            f"- peak 위치 변화: {1000*v['max_peak_location_change_m']:.6f} mm (기준 2 mm)",
            f"- CPU/GPU peak 상대차: {row['cpu_gpu_peak_relative_difference']:.3g}",
            f"- CPU/GPU q 최대 절대차: {row['cpu_gpu']['max_absolute_q']:.3g}",
            f"- 화면 형상/정밀 벽 최대 거리: {1e6*row['profile_wall_hausdorff_m']:.3f} μm (기준 20 μm)",
            '', '| 보간 | 경로 수 | 기준 peak | 선택 peak |','|---|---:|---:|---:|']
        for level in row['levels']:lines.append(f"| {level['method']} | {level['tubes']} | {level['baseline_peak']:.9f} | {level['selected_peak']:.9f} |")
        lines+=['']
    lines+=['## 보존과 표시','',
        '원래 campaign, 기준 입력, 원래 manifest/selected_refined, viewer data.json·legs.json·기존 단계 프로필은 변경하지 않았다. '
        '새 결과는 별도 study의 QualifiedStore에 체크섬과 함께 영구 보존했고, 뷰어에는 추가 정밀검증 결과라는 별도 표시로 연결했다. '
        '기존 단계들은 여전히 320 경로 탐색 기록이며, 추가 최종 결과는 1280 경로이다.','',
        '## 해석 범위','',
        '선택 후보는 이번에 평가한 보관 후보 목록 내에서 모든 검증을 통과한 최저 1280 peak 후보이다. '
        '근접한 후보 사이의 차이는 격자·보간 오차보다 작을 수 있으므로 이 순위를 물리적으로 유의한 우열로 단정하지 않는다. '
        '연속 형상 변수 공간의 전역 최적성 증명은 아니다. 전력 보존·수치 수렴은 사용한 축약 수송 모델의 검증이며 '
        '실험값 일치나 완전한 플라스마 물리 검증을 뜻하지 않는다. OMP 100 mm는 이번 추가 검증 대상이 아니다.','',
        '![형상과 열유속 비교](omp-qualified-comparison.png)','',
        f"Study fingerprint: `{report['config']['fingerprint']}`",'',
        f"추가 검증 실행 시간: {report['elapsed_seconds']:.1f} s (초기 입력 해시 점검 제외).",'',
        '기계 판독 검수: viewer/qualified-final-verification.json. 원시 순위 및 전체 검증: study/summary.json 및 R020/R050. GitHub에는 push하지 않았다.']
    (dest/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    save(dest/'report_sha256.json',{p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name!='report_sha256.json'})
    print(dest/'REPORT.md')


if __name__=='__main__':main()
