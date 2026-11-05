"""Reproducible figures and Korean summary from qualified campaign evidence."""
import argparse,json,csv,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from d3d_optimization.broadening_study import load_json,digest
from d3d_optimization.campaign_state import public_save as save
from d3d_optimization.campaign_incumbents import QualifiedStore,context


def selection(folder,expected):
    """Recover the canonical bundle even after a convenience-pointer crash."""
    folder=Path(folder);record=QualifiedStore(folder,expected).best()
    if record is not None:return record,folder/record['artifact_bundle']/'selected_refined.npz'
    path=folder/'manifest.json';m=load_json(path) if path.exists() else {}
    if m.get('artifact_bundle'):raise ValueError('qualified manifest references a missing bundle')
    return m,folder/'selected_refined.npz'  # Read-only legacy report support.


def run(folder):
    folder=Path(folder);config=load_json(folder/'campaign.json');rows=[];manifests={};sources={};selected_paths={};evidence_hashes={}
    for source in ['IMP','OMP']:
        path=folder/source/'worker_status.json';sources[source]=load_json(path) if path.exists() else {'status':'not_started'}
        hard_stop=folder/source/'hard_deadline.json'
        if hard_stop.exists():sources[source]=dict(sources[source],status='deadline_stopped',hard_deadline=load_json(hard_stop))
        for mm in [20,50,100]:
            path=folder/source/f'R{mm:03d}'/'manifest.json';progress=folder/source/f'R{mm:03d}'/'progress.json'
            m,selected_paths[(source,mm)]=selection(path.parent,context(config,source,mm))
            p=load_json(progress) if progress.exists() else {}
            qualified=m.get('status')=='qualified' and m.get('validation',{}).get('verdict',{}).get('resolution_screen_passed')
            status='qualified' if qualified else m.get('status','searching')
            if not qualified and p.get('best') and sources[source]['status'] in ('searching','validating','starting'):
                status='validation_pending'
            if not qualified and p.get('best') and sources[source]['status']=='deadline_stopped' and len(m.get('qualification_attempts',[]))<4:
                status='validation_incomplete_deadline'
            row=dict(source=source,radius_mm=mm,status=status,
                baseline_peak_MW_m2=m.get('baseline',{}).get('peak_MW_m2') if qualified else None,
                selected_peak_MW_m2=m.get('selected',{}).get('peak_MW_m2') if qualified else None,
                reduction_percent=100*(1-m['selected']['peak_ratio']) if qualified else None,
                proposals=p.get('total_proposals',0),physics_evaluations=p.get('physics_evaluations',0),
                restarts=p.get('restart',0),coarse_best_ratio=(p.get('best') or {}).get('objective'),
                global_optimum_proven=False)
            rows.append(row)
            if qualified:
                manifests[(source,mm)]=m
                evidence=(path.parent/m['artifact_bundle']/'manifest.json') if m.get('artifact_bundle') else path
                evidence_hashes[evidence.relative_to(folder).as_posix()]=digest(evidence)
    save(folder/'summary.json',dict(rows=rows,worker_status=sources,epoch=time.time(),fingerprint=config['fingerprint'],
        source_manifest_hashes=evidence_hashes,
        author_script_sha256=digest(Path(__file__)),global_optimum_proven=False))
    with (folder/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors={20:'#176B9C',50:'#D68A16',100:'#9B426B'}
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,source in zip(axes,['IMP','OMP']):
        for mm in [20,50,100]:
            path=folder/source/f'R{mm:03d}'/'evaluations.jsonl';x=[];y=[];best=np.inf
            if path.exists():
                with path.open(encoding='utf-8') as f:
                    for line in f:
                        try:r=json.loads(line)
                        except ValueError:continue  # interrupted last log line is not a checkpoint
                        if r['feasible'] and r['objective']<best:
                            best=r['objective'];x.append((r['epoch']-config['created_epoch'])/3600);y.append(best)
            if x:ax.step(x,y,where='post',label=f'{mm} mm',color=colors[mm])
        ax.axhline(1,color='#333333',linestyle='--',linewidth=1,label='Baseline')
        ax.set(title=source,xlabel='Elapsed wall time [h]',ylabel='Best search peak / baseline (320 tubes)')
        ax.grid(alpha=.15);ax.legend()
    power_mw=config['reference']['power_W']/1e6
    fig.suptitle(f'Independent {power_mw:g} MW sources — search progress (final qualification is separate)')
    fig.savefig(folder/'search_progress.png',dpi=170);plt.close(fig)
    if manifests:
        maximum=max(max(m['baseline']['peak_MW_m2'],m['selected']['peak_MW_m2']) for m in manifests.values())
        fig,axes=plt.subplots(2,4,figsize=(17,8),layout='constrained');norm=Normalize(0,np.ceil(maximum*10)/10)
        for row_index,source in enumerate(['IMP','OMP']):
            for col,mm in enumerate([None,20,50,100]):
                ax=axes[row_index,col];label=f'{source} baseline' if mm is None else f'{source}, Rmin {mm} mm'
                path=folder/source/'baseline_refined.npz' if mm is None else selected_paths[(source,mm)]
                certified=mm is None or (source,mm) in manifests
                if path.exists() and certified:
                    with np.load(path) as d:
                        wall=d['receiving_wall'];power=d['cell_power']*config['reference']['power_W']/1e6
                        q=power/np.diff(d['area_edges']);segments=np.stack([wall[:-1],wall[1:]],axis=1)
                        ax.add_collection(LineCollection(segments,array=q,cmap='magma',norm=norm,linewidths=2.5))
                    m=next((m for (s,_),m in manifests.items() if s==source),None) if mm is None else manifests[(source,mm)]
                    peak=m['baseline' if mm is None else 'selected']['peak_MW_m2'] if m else None
                    if peak is not None:label+=f'\nqmax = {peak:.4f} MW/m²'
                else:
                    row=next((r for r in rows if r['source']==source and r['radius_mm']==mm),{})
                    label_status={'validation_pending':'Qualification pending',
                                  'validation_incomplete_deadline':'Validation budget ended'}.get(row.get('status'),'No qualified candidate')
                    ax.text(.5,.5,label_status,ha='center',va='center',transform=ax.transAxes)
                ax.set(title=label,xlabel='R [m]',ylabel='Z [m]')
                ax.set_xlim((.99,1.19) if source=='IMP' else (1.12,1.72));ax.set_ylim((-1.37,-1.10) if source=='IMP' else (-1.46,-1.20))
                ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='magma'),ax=axes.ravel().tolist(),shrink=.8,label='Cell-average heat flux [MW/m²]')
        fig.suptitle('Qualified geometry at 1280 tubes — each source optimized separately')
        fig.savefig(folder/'qualified_geometry.png',dpi=170);plt.close(fig)
    deadline=datetime.fromtimestamp(config['deadline_epoch'],timezone.utc).astimezone().isoformat(timespec='seconds')
    lines=['# IMP·OMP 접선·곡률 제약 전역 탐색', '',
        f'- 종료 시한: {deadline}', f'- 열원별로 {power_mw:g} MW를 부여한 독립 실험입니다.',
        '- 기준 D3D, 입력 전력과 수송 조건은 고정했습니다. 표의 최종 peak는 상류 flux tube 1280개 결과입니다.',
        '- 전역 최적성 증명은 아닙니다. 정해진 곡선 계열과 계수 범위에서 찾은 최선 후보입니다.',
        '- IMP 공통 접선은 두 구간에 대해 정규화 내적 ≥0.2인 방향 범위에서 탐색합니다.',
        '- IMP의 3개 연결, OMP 바닥 왼쪽 연결에 G1 조건을 적용합니다. 곡률 제한은 양쪽 부호를 허용합니다.',
        '- OMP 오른쪽 개구부는 기존 최소 간격을 유지합니다. 원래 고정된 다른 공학 코너는 그대로입니다.', '',
        '| 영역 | Rmin (mm) | 상태 | 기준 peak | 선택 peak (MW/m²) | 감소율 | 제안 / 열계산 | 재시작 |',
        '|---|---:|---|---:|---:|---:|---:|---:|']
    for r in rows:
        fmt=lambda x:'—' if x is None else f'{x:.4f}'
        lines.append(f"| {r['source']} | {r['radius_mm']} | {r['status']} | {fmt(r['baseline_peak_MW_m2'])} | {fmt(r['selected_peak_MW_m2'])} | {fmt(r['reduction_percent'])}% | {r['proposals']} / {r['physics_evaluations']} | {r['restarts']} |")
    lines+=['','후보 미발견은 불가능성 증명이 아닙니다. 조건을 완화하지 않고 제약 위반 기록을 보존했습니다.',
        '정밀 검증은 전력 보존·전체 첫 충돌·독립 열적분·도착 위치·열유속 CDF·격자 수렴과 GPU/CPU 분포 일치를 포함합니다.',
        '수송은 기존의 단방향 열수송 및 전력 보존형 벽 broadening 축약 모델입니다. 측정값 일치는 이 연구의 판정 조건이 아닙니다.',
        '', '![탐색 이력](search_progress.png)']
    if manifests:lines+=['','![검증된 형상](qualified_geometry.png)']
    (folder/'REPORT.ko.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(folder),qualified=len(manifests),rows=rows),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
