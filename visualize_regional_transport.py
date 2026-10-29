"""Scientific figures and data exports for independent regional searches."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from d3d_optimization.broadening_study import load_json,save,digest

COLORS={'baseline':'#2876b5','selected':'#d66a27'}


def profile(a,power):
    wall=a['receiving_wall'];area=np.diff(a['area_edges']);ds=np.linalg.norm(np.diff(wall,axis=0),axis=1)
    arc=np.r_[0.,np.cumsum(ds)];s=(arc[:-1]+arc[1:])/2
    return wall,s,a['cell_power']/area*power/1e6


def render(output):
    output=Path(output);loaded={};summary=[]
    for source in ['IMP','OMP']:
        folder=output/'qualified'/source if (output/'qualified'/source/'manifest.json').exists() else output/source
        relative=folder.relative_to(output).as_posix()
        m=load_json(folder/'manifest.json');v=load_json(folder/'validation.json');rv=load_json(folder/'regional_validation.json')
        if not v['complete'] or not rv['passed']:raise ValueError('complete validation and regional invariance required')
        rows=v['levels'][2]['rows'];arrays={n:np.load(folder/(n+'_refined.npz')) for n in ['baseline','selected']}
        history=[json.loads(line) for line in (folder/'search/evaluations.jsonl').read_text().splitlines()]
        # Flat inspectable history, including rejected candidates and exact regional coefficients.
        header=['case','stage','feasible','peak_MW_m2','peak_ratio','elapsed_seconds','opposite_region_displacement_m','errors']
        coeff_names=[f'{name}_c{k+1}_m' for name in m['regional_restriction'] for k in range(4)]
        with (folder/'optimization_history.csv').open('w',newline='',encoding='utf-8') as f:
            writer=csv.writer(f);writer.writerow(header+coeff_names)
            for r in history:writer.writerow([json.dumps(r[k]) if k=='errors' else r.get(k,'') for k in header]+r['coefficients_m'])
        m['delivery_relative_folder']=relative
        loaded[source]=(m,v,rows,arrays,history)
        p=m['reference']['power_W'];b,s=rows['baseline'],rows['selected']
        summary.append(dict(source=source,baseline_peak_MW_m2=b['peak_MW_m2'],selected_peak_MW_m2=s['peak_MW_m2'],
            reduction_percent=100*(1-s['peak_ratio']),confirmed=v['verdict']['improvement_confirmed'],
            baseline_effective_area_m2=b['ledger']['deposited']*p/(b['peak_MW_m2']*1e6),
            selected_effective_area_m2=s['ledger']['deposited']*p/(s['peak_MW_m2']*1e6),
            baseline_L_parallel_mean_m=b['thermal']['length_m']['power_weighted_mean'],
            selected_L_parallel_mean_m=s['thermal']['length_m']['power_weighted_mean'],
            baseline_sigma_mean_mm=b['sigma_physical_mean_m']*1000,selected_sigma_mean_mm=s['sigma_physical_mean_m']*1000,
            max_displacement_mm=s['max_wall_displacement_m']*1000,opposite_region_error_m=rv['opposite_and_fixed_wall_max_error_m'],
            selected_deposited_power_W=s['ledger']['deposited']*p,
            selected_unhit_power_W=s['ledger']['unhit']*p,
            selected_unresolved_power_W=s['ledger']['unresolved']*p,
            selected_power_balance_error_fraction=s['ledger']['balance_error'],
            max_peak_ratio_change=v['verdict'].get('max_peak_ratio_change'),
            evaluations=m['search']['evaluated'],feasible_evaluations=m['search']['feasible_evaluations'],
            optimization_seconds=m['elapsed_seconds'],validation_level_seconds_excluding_wall_only=sum(l['elapsed_seconds'] for l in v['levels']),backend=m['backend'],
            result_folder=relative,qualified_alternative='qualification' in m))
        fig,ax=plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for name in ['baseline','selected']:
            a=arrays[name];color=COLORS[name];style='--' if name=='baseline' else '-'
            ax[0,0].plot(a['wall'][:,0],a['wall'][:,1],style,color=color,label=name,lw=1.8)
            w,x,q=profile(a,p);ax[0,1].plot(x,q,style,color=color,label=name,lw=1.6)
            significant=x[q>.01*b['peak_MW_m2']]
            if name=='baseline':extent=[significant.min()-.03,significant.max()+.03]
            else:extent=[min(extent[0],significant.min()-.03),max(extent[1],significant.max()+.03)]
            ax[1,1].plot(a['source_x_m']*1000,a['sigma_physical_m']*1000,style,color=color,label=name)
        bounds=((.99,1.19),(-1.38,-1.12)) if source=='IMP' else ((1.14,1.71),(-1.41,-1.18))
        ax[0,0].set(xlim=bounds[0],ylim=bounds[1],xlabel='R [m]',ylabel='Z [m]',title=f'{source} design region; other region fixed')
        ax[0,0].set_aspect('equal');ax[0,0].legend()
        ax[0,1].set(xlim=extent,xlabel='Connected wall arclength [m]',ylabel='Cell-average q [MW/m²]',title='Arriving heat-flux distribution');ax[0,1].legend()
        cases=np.array([r['case'] for r in history]);best=[];incumbent=1.
        for r in history:
            if r['feasible']:incumbent=min(incumbent,r.get('peak_ratio',1.))
            best.append(incumbent)
        good=[r for r in history if r['feasible'] and 'peak_ratio' in r]
        ax[1,0].scatter([r['case'] for r in good],[r['peak_ratio'] for r in good],s=9,color='#b5bbc2',label='feasible proposal')
        ax[1,0].step(cases,best,where='post',color=COLORS['selected'],lw=2,label='original search minimum')
        if 'qualification' in m:
            chosen=next(r for r in history if r['case']==m['qualification']['case'])
            ax[1,0].scatter([chosen['case']],[chosen['peak_ratio']],s=95,marker='*',color='#222222',label='selection after validation',zorder=5)
        ax[1,0].axhline(1,color='#444444',ls='--',lw=1);ax[1,0].set(xlabel='Unique evaluated candidate',ylabel='Peak / baseline peak',title='Optimization history (search resolution)');ax[1,0].legend()
        vals=np.concatenate([a['sigma_physical_m'][a['source_x_m']<=.02]*1000 for a in arrays.values()]);vals=vals[np.isfinite(vals)&(vals>0)]
        ax[1,1].set(xlim=(0,20),ylim=(vals.min()*.8,vals.max()*1.25),yscale='log',xlabel=f'{source} source offset [mm]',ylabel='Effective footprint σ [mm]',title='Transport width at identical inputs');ax[1,1].legend()
        for aa in ax.ravel():aa.grid(alpha=.2)
        fig.suptitle(f"{source}: {b['peak_MW_m2']:.4f} → {s['peak_MW_m2']:.4f} MW/m² | reduction {100*(1-s['peak_ratio']):.2f}%",fontsize=15)
        fig.savefig(folder/'optimization_summary.png',dpi=180);plt.close(fig)
        # Refined cell/donor tables are separate from search-resolution exports.
        for n,a in arrays.items():
            w,x,q=profile(a,p);mid=(w[:-1]+w[1:])/2
            np.savetxt(folder/(n+'_profile_refined.csv'),np.c_[x,mid,q,a['cell_power']*p,np.diff(a['area_edges'])],delimiter=',',
                header='wall_arclength_m,R_m,Z_m,cell_average_q_MW_m2,cell_power_W,cell_area_m2',comments='')
            np.savetxt(folder/(n+'_donors_refined.csv'),np.c_[a['source_x_m'],a['hit_points'],a['thermal_length_m'],a['thermal_Tu_eV'],
                a['thermal_Tt_eV'],a['sigma_physical_m'],a['donor_region_power'].sum(axis=1)*p],delimiter=',',
                header='offset_m,hit_R_m,hit_Z_m,L_parallel_m,Tu_eV,Tt_eV,sigma_physical_m,delivered_donor_power_W',comments='')
    # Same color scale for all four panels permits an absolute peak comparison.
    vmax=max(row['baseline_peak_MW_m2'] for row in summary)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for i,source in enumerate(['IMP','OMP']):
        m,v,rows,arrays,history=loaded[source]
        for j,name in enumerate(['baseline','selected']):
            ax=axes[i,j];a=arrays[name];w,x,q=profile(a,m['reference']['power_W'])
            lc=LineCollection(np.stack([w[:-1],w[1:]],axis=1),cmap='inferno',norm=plt.Normalize(0,vmax),linewidth=4)
            lc.set_array(q);ax.add_collection(lc)
            ax.set(xlim=(.98,1.72),ylim=(-1.43,-1.10),xlabel='R [m]',ylabel='Z [m]',
                title=f"{source} — {name} | qmax {rows[name]['peak_MW_m2']:.4f} MW/m²")
            ax.set_aspect('equal');ax.grid(alpha=.15)
    fig.colorbar(lc,ax=axes,label='Cell-average heat flux [MW/m²]',shrink=.85)
    fig.suptitle('Independent IMP and OMP experiments — 1 MW assigned to each source separately',fontsize=14)
    fig.savefig(output/'heat_flux_comparison.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,source in enumerate(['IMP','OMP']):
        for name,a in loaded[source][3].items():
            for j,(key,scale,label) in enumerate([
                ('thermal_length_m',1.,'Parallel connection length [m]'),
                ('sigma_physical_m',1000.,'Effective footprint σ [mm]')]):
                axes[i,j].plot(a['source_x_m']*1000,a[key]*scale,
                    '--' if name=='baseline' else '-',color=COLORS[name],label=name)
                axes[i,j].set(xlim=(0,20),xlabel=f'{source} source offset [mm]',ylabel=label,
                    title=f'{source}: '+('path including toroidal winding' if j==0 else 'transport footprint'))
                if j==1:axes[i,j].set_yscale('log')
                axes[i,j].grid(alpha=.2);axes[i,j].legend()
        visible=np.concatenate([a['sigma_physical_m'][a['source_x_m']<=.02]*1000
            for a in loaded[source][3].values()])
        visible=visible[np.isfinite(visible)&(visible>0)]
        axes[i,1].set_ylim(visible.min()*.8,visible.max()*1.25)
    fig.suptitle('Path length and broadening at unchanged operating inputs',fontsize=14)
    fig.savefig(output/'transport_comparison.png',dpi=180);plt.close(fig)
    with (output/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    gpu=load_json(output/'gpu_preflight.json')
    lines=['# IMP·OMP 영역별 독립 형상 최적화','',
        '기준 D3D 형상과 자기 평형을 보존하고, 같은 입력 전력에서 도착 열유속을 넓게 분산해 peak를 낮추는 형상을 탐색했습니다.',
        'IMP 열원에서는 imp_sloped·imp_vertical만, OMP 열원에서는 omp_entry·shelf·omp_connector·floor만 변형했습니다. 각 실행의 반대 영역은 고정입니다.',
        '**두 열원 각각 1 MW를 부여한 독립 실험입니다. 두 열원을 동시에 켠 합성 형상의 최적화 결과는 아닙니다.**','',
        '입력: 각 소스 P=1 MW, n=1e19 m⁻³, χ⊥=0.3 m²/s, λq=3 mm, γ=7, κ₀=2000. 소스의 유한 계산 구간은 각 manifest에 기록했습니다.','',
        '| 소스 | 기준 peak | 선택 peak | 감소율 | 수치 검증 | 고유 후보 | 탐색 시간 |','|---|---:|---:|---:|---|---:|---:|']
    for r in summary:lines.append(f"| {r['source']} | {r['baseline_peak_MW_m2']:.6f} | {r['selected_peak_MW_m2']:.6f} | {r['reduction_percent']:.3f}% | {r['confirmed']} | {r['evaluations']} | {r['optimization_seconds']:.1f} s |")
    lines.extend(['','열유속 단위는 MW/m²입니다. 최종 수치는 cubic/1280 donor 결과이며, 320/640/1280 및 quintic640 보간·격자 검증을 통과한 경우에만 개선 확인을 표시합니다. 시간은 두 탐색을 병렬 실행했을 때 각 프로세스의 경과시간으로, 초기 계산 준비를 포함하고 후속 후보 선별과 정밀 검증 시간은 제외합니다.','',
        f"GPU: RTX 3090의 FP64 Gaussian 합산 사용. 같은 OMP320 후보 세 개의 중앙 평가시간 비교에서 약 {gpu['gpu_speedup']:.2f}배 빨랐습니다. CPU/GPU 전 분포 차이는 약 5e-16(peak 대비), 최대값 차이는 0이었습니다. 이 가속비는 전체 최적화의 가속비를 뜻하지 않습니다.",
        'CDF 전력·형상 교차·열 상태·연속 최대점 정밀 탐색 및 최종 검증은 CPU로 수행했습니다. GPU 패키지, 캐시, 연구 데이터와 코드/입력 스냅샷은 D 드라이브에 저장했습니다.','',
        '## 형상과 분포','', '![열유속 지도](heat_flux_comparison.png)',
        '', '![연결 길이와 수송 폭](transport_comparison.png)'])
    for r in summary:
        source=r['source'];m,v,rows,arrays,h=loaded[source]
        rel=m['delivery_relative_folder']
        lines.extend(['',f'### {source}','',f"![{source} 형상·분포·탐색 이력]({rel}/optimization_summary.png)",'',
          f"반대 영역/고정 벽의 최대 좌표 변화: {r['opposite_region_error_m']:.1e} m. 최대 설계 변위: {r['max_displacement_mm']:.3f} mm.",
          f"도착 전력: {r['selected_deposited_power_W']/1e6:.6f} MW; 미도착/미해결 전력: {r['selected_unhit_power_W']:.3g}/{r['selected_unresolved_power_W']:.3g} W; 전력 장부 상대 오차: {r['selected_power_balance_error_fraction']:.2e}.",
          f"유효 가열 면적 P도착/qmax: {r['baseline_effective_area_m2']:.4f} → {r['selected_effective_area_m2']:.4f} m². 이는 FWHM이나 특정 열유속 문턱의 면적과 다른 지표입니다.",
          f"전력가중 평균 연결 길이: {r['baseline_L_parallel_mean_m']:.3f} → {r['selected_L_parallel_mean_m']:.3f} m; 유효 footprint 폭: {r['baseline_sigma_mean_mm']:.3f} → {r['selected_sigma_mean_mm']:.3f} mm.",
          f"최대 peak 상대변화(격자/보간): {100*v['verdict'].get('max_absolute_peak_relative_change',float('nan')):.4f}%; 최대점 위치 변화: {1000*v['verdict'].get('max_peak_location_change_m',float('nan')):.4f} mm.",
          f"격자/보간 간 peak 비율 차이: {100*v['verdict'].get('max_peak_ratio_change',float('nan')):.4f}%p (허용치 0.2%p).",
          f"[선택 전체 벽 좌표]({rel}/selected_wall_refined.csv) · [선택 도착 분포]({rel}/selected_profile_refined.csv) · [후보별 탐색 데이터]({rel}/optimization_history.csv) · [입력·최종 계수]({rel}/manifest.json)"])
        if 'qualification' in m:
            screen=load_json(Path(m['qualification']['screen_path']))
            original=load_json(Path(m['original_search_best']['validation_path']))
            raw_peak=original['levels'][2]['rows']['selected']['peak_MW_m2']
            lines.extend(['',f"원래 탐색 최저 peak 후보의 정밀 격자값은 {raw_peak:.6f} MW/m²였으나, 보간법에 따른 peak 비율 차이가 허용치 0.2%p를 넘어 개선 확정 판정을 받지 못했습니다. 원래 결과는 `{source}/`에 보존했습니다. 기준을 바꾸지 않고 탐색 로그에서 순위에 따라 뽑은 {len(screen['candidates'])}개 후보를 검사하여 그중 통과 후보의 peak가 가장 낮은 case {m['qualification']['case']}(원래 탐색 순위 {m['qualification']['rank']})를 선택했으며, 이 후보에 전체 검증을 다시 수행했습니다. 탐색 이력의 최저선과 검증 후 선택점이 다른 이유입니다."])
    lines.extend(['','## 데이터 사용과 한계','',
        '- 각 영역 폴더의 `selected_<표면명>_refined.csv`가 해당 표면의 최종 R,Z 좌표입니다. `selected_wall_refined.csv`는 반대 영역을 고정한 전체 벽입니다.',
        '- `*_refined.npz`에는 전 도착 열유속 샘플, 셀 전력, 자기장 기울기, 연결 길이, 열 상태와 donor→표면 전달 배열이 있습니다. CSV는 검사·후처리용입니다.',
        '- 입력 및 수송 가정은 기존과 동일합니다. 주기 Gaussian을 연결된 벽의 면적 좌표에 분배하는 대리모델이며, 실제 2D 수송 방정식이나 방사 손실을 풀지 않습니다.',
        '- 토로이달 감김을 포함한 연결 길이를 적분하지만, 이번 감소율을 연결 길이 단독 효과로 분해한 결과는 아닙니다.',
        '- 목적함수는 연결 길이나 평균 σ 자체를 키우는 것이 아니라 도착 peak를 낮추는 것입니다. 입사 위치, 표면 기하, 면적 및 수송 폭이 함께 바뀌므로 평균 연결 길이 또는 평균 σ가 감소한 후보도 peak를 낮출 수 있습니다.',
        '- 분포 그림의 가로축은 각 형상을 따른 벽 호길이입니다. 같은 호길이가 두 형상에서 정확히 같은 R,Z 위치를 뜻하지는 않습니다. 최대값은 연속 분포에서 정밀 탐색하며, 열유속 지도와 곡선은 셀 평균값입니다.',
        '- 제한된 탐색 예산에서 찾은 적합 후보입니다. 전역 최적해 또는 두 독립 최적 형상의 조합이 동시 두 열원에 최적이라는 주장은 하지 않습니다.',
        '- 두 결과를 한 형상으로 합치려면 합성 벽에서 두 열원의 열유속을 함께 다시 평가해야 합니다. 이번 데이터에는 각 영역의 독립 결과를 분리해 보존했습니다.'])
    (output/'REPORT.ko.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    save(output/'summary.json',dict(results=summary,gpu_preflight=gpu,plot_script_sha256=digest(Path(__file__))))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();print(render(a.output))
