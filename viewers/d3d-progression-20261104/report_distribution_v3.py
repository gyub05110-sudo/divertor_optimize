"""Read-only final-study verification and standalone report/plots."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qualify_omp_stages import read,sha,save,check_preservation
from publish_distribution_v3 import OUT,HERE,STUDY
from distribution_policy_v3_final import from_npz

def main():
    marker=read(OUT/'COMPLETED.json');report=read(OUT/'summary.json');config=read(OUT/'config.json')
    assert marker['summary_sha256']==sha(OUT/'summary.json') and marker['artifact_index_sha256']==sha(OUT/'artifact_sha256.json')
    for rel,h in read(OUT/'artifact_sha256.json').items():assert sha(OUT/rel)==h,rel
    assert check_preservation(config['protected_sha256'])['unchanged'] and check_preservation(config['code_sha256'])['unchanged']
    rows=report['candidates'];assert len(rows)==20 and sum(r['radius_mm']==20 for r in rows)==19
    assert len({r['identity'] for r in rows})==20
    for mm in (20,50):
        passed=sorted([r for r in rows if r['radius_mm']==mm and r['passed']],key=lambda r:(r['peak'],r['identity']))
        assert report['results'][str(mm)]['best']==(passed[0] if passed else None)
    out=OUT/'report';out.mkdir(exist_ok=True)
    lines=['# 분포 수렴 v3 재검증 결과','',
        '고정 D3D 기준, 독립 1 MW, 기존 수송 입력을 유지했다. 새 형상 탐색 없이 보관 후보를 재평가했다. 원래 연구와 v2 검증/최종 형상/뷰어 데이터는 변경하지 않았다.','',
        '## 판정 변경','',
        '- 최대점 위치 2 mm는 합격조건에서 제외하고 진단값으로 보존한다.',
        '- 동일한 물리적 벽 좌표에서 저장된 연속 q 샘플의 구간별 선형 재구성을 비교한다. 봉우리 정렬/분포 이동/전력 재정규화는 하지 않는다.',
        '- 필수: 최대 위치별 차이 / 비교 peak ≤1%, 면적 가중 L1 / 평균 전력 ≤0.5%, A2 변화 ≤0.5%, 전력 변화 ≤0.5%, 벽 대응/길이 차이 ≤50 μm.',
        '- 실제 셀 전력 합 오차 ≤2e-12와 재구성 전력 오차 ≤0.1%를 별도로 검사한다. A2=P²/∫q²dA이며 dA=2πR ds. 원래 셀 평균 A2도 수렴해야 한다.',
        '- 기존 peak·ratio·영역별 전력·CDF·개선 여유·독립 first-hit/열적/형상·벽 세분화 조건은 유지했다. 기준과 후보 모두 fine 3쌍 및 벽 전용 1쌍을 검사했다.',
        '- 기준값은 연구용 수치 허용오차이며 제작/물리 표준이 아니다. 후보를 통과시키려고 사후 완화하지 않았다.','',
        '## 평가 결과','',
        '| Rmin | 원래 단계 | peak MW/m² | 이전 기준 | v3 | 최대 q 차이 / peak | L1 | A2 변화 |',
        '|---:|---:|---:|---|---|---:|---:|---:|']
    for r in sorted(rows,key=lambda x:(x['radius_mm'],x['stage'])):
        lines.append(f"| {r['radius_mm']} | {r['stage']} | {r['peak']:.9f} | {'통과' if r['previous_pass'] else '미통과'} | {'통과' if r['passed'] else '미통과'} | {100*r['max_linf_relative']:.4f}% | {100*r['max_l1_relative']:.4f}% | {100*r['max_a2_relative']:.4f}% |")
    lines+=['','OMP20은 기존 정밀평가 16개(탈락15개+이전 최종1개)에 초기3개를 추가해 보관 후보19개 전체를 검사했다. OMP50은 53개 원래 fine 순위 중 해당 반경 최저 peak 후보가 통과하므로, 더 높은33개 후보의 전체 검증을 반복하지 않았다.','',
        '## 적대적 검증과 한계','',
        '비슷한 쌍봉우리의 높이 순위 교체와 동일한 매끄러운 함수를 다른 셀 폭으로 평균한 경우는 통과한다. 실제 분포 이동, 봉우리 소실, 전역 최대값보다 낮은 좁은 부봉우리, 전력 손실, 벽 좌표 반전/변경, 독립 검사 실패는 거부한다. 절댓값 오차의 영점에서 적분 구간을 분할한다.',
        '원래 서로 다른 셀 평균을 계단함수로 직접 비교하면 평균 구간 차이가 Linf 오차에 섞인다. 이 프로토타입과 최초 검사 기록도 보존하고, 최종 판정에서는 같은 위치의 연속 샘플 재구성을 사용했다. 허용오차는 동일하다.',
        '저장 샘플 사이의 완전한 연속 함수 오차를 수학적으로 보증하지 않는다. 독립 dense peak 검사는 유지하지만 표본에서 누락된 임의의 국소 구조까지 배제하는 인증은 아니다. 보간/격자 변화에 대한 수치 일관성이며 실험 일치나 전역 최적성 증명이 아니다.',
        'Linf만 절반/두 배로 바꾼 민감도(나머지 검사 통과를 뜻하지 않음): '+json.dumps(report['sensitivity']),
        '', '기준 형상의 새 분포 검사: '+str(report['baseline']['passed']),
        '', 'GitHub push 없음. 새 검증/출판 시간: '+f"{report['elapsed_seconds']:.1f} s (초기3개 추가 열수송 검증 시간 별도).",'']
    # Resolve the actual largest pointwise difference of the former raw-best candidate.
    identity='0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1'
    p=STUDY/f'R020/validation/{identity}';a=from_npz(p/'cubic_640.npz');b=from_npz(p/'quintic_640.npz')
    ua,qa=a['continuous'];ub,qb=b['continuous'];u=np.unique(np.r_[ua,ub]);va=np.interp(u,ua,qa);vb=np.interp(u,ub,qb)
    j=int(np.argmax(abs(va-vb)));R=np.interp(u,a['u'],a['wall'][:,0]);Z=np.interp(u,a['u'],a['wall'][:,1]);arc=u*a['length']
    diagnostic=dict(identity=identity,R_m=float(R[j]),Z_m=float(Z[j]),q_cubic640=float(va[j]),q_quintic640=float(vb[j]),absolute_difference=float(abs(va[j]-vb[j])))
    save(out/'largest-local-difference.json',diagnostic)
    fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained',sharex=True)
    mask=(abs(arc-arc[j])<.02)
    axes[0].plot(arc[mask]-arc[j],va[mask],label='cubic 640');axes[0].plot(arc[mask]-arc[j],vb[mask],label='quintic 640');axes[0].set(ylabel='q [MW/m²]',title='OMP20 stage18: local distribution difference beyond argmax switching');axes[0].legend()
    axes[1].plot(arc[mask]-arc[j],100*abs(va[mask]-vb[mask])/max(qa.max(),qb.max()));axes[1].axhline(1,color='red',ls='--',label='v3 tolerance: 1%');axes[1].set(xlabel='Wall arc distance from largest discrepancy [m]',ylabel='|difference| / global peak [%]');axes[1].legend()
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(out/'local-distribution-difference.png',dpi=180);plt.close(fig)
    lines+=['## 최대 국소 차이', '',json.dumps(diagnostic,ensure_ascii=False),'','![국소 q 차이](local-distribution-difference.png)','',
            '최종 결과 인덱스: ../viewer/distribution-finals.json. 이전 기준20 mm 최종 형상은 파일을 보존하되 v3 미통과로 구분한다.']
    (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    save(out/'checksums.json',{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='checksums.json'})
    save(HERE/'distribution-study-verification.json',dict(passed=True,all_core_hashes_verified=True,protected_inputs_unchanged=True,case_count=len(rows),results=report['results']))
    print('VERIFIED',report['results'])

if __name__=='__main__':main()
