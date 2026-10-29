"""Compare verified baseline/candidate heat samples, without running a search."""
import csv
import json
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from divertor_optimization.d3d_multisurface import MultiSurfaceEvaluator
from divertor_optimization.d3d_frozen_case import load_case, separatrix, pretrace, sha256, WIDTHS_MM

p=Path(__file__).resolve().parent
out=p/'comparison';out.mkdir(exist_ok=True)
m=json.loads((p/'manifest.json').read_text(encoding='utf-8'))
assert m['numerically_verified_feasible']
assert sha256(m['shot_path'])==m['shot_sha256']
assert sha256(p/'DIIID_geom.json')==m['config']['geometry_sha256']
g=json.loads((p/'DIIID_geom.json').read_text())
wall=np.asarray(g['limiter']);wall=np.vstack([wall,wall[0]])
field,lcfs,axis=load_case(m['shot_path'],m['time_ms'],'cubic')
r_sep,_,_=separatrix(field,lcfs,axis)
edges,paths,inside,_,_=pretrace(field,r_sep,axis[1],wall,320,.002)
e=MultiSurfaceEvaluator(field,lcfs,wall,paths,inside,edges,r_sep,axis[1],m['config']['surfaces'],
                        vessel_solids=shapely.union_all([Polygon(v['pts']) for v in g['vv']]))
candidate=e.evaluate(m['coefficients_m'])
assert candidate['feasible']
results=[e.baseline,candidate]
labels=['기준 형상','최적화 형상']
colors=['#2166ac','#e07a1f']
peaks=np.array([[max(h.q) for h in r['maps']] for r in results])
reviewed=list(csv.DictReader((p/'width_scan.csv').open()))
np.testing.assert_allclose(peaks.ravel(),[float(r['all_wall_peak_q_over_P']) for r in reviewed],rtol=1e-10)
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False,'font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})

samples=[]
for label,result in zip(labels,results):
    for mm,h in zip(WIDTHS_MM,result['maps']):
        for i,seg in enumerate(h.segment):
            component='미포착' if seg<0 else '고정 표면'
            for name,ids in result['surface_segments'].items():
                if seg in ids:component=name
            samples.append({'geometry':label,'lambda_q_mm':int(mm),'tube':i,'component':component,
                            'R_m':h.points[i,0],'Z_m':h.points[i,1],'q_over_P_outer_per_m2':h.q[i]})
with (out/'heat_samples.csv').open('w',newline='',encoding='utf-8-sig') as f:
    writer=csv.DictWriter(f,fieldnames=list(samples[0]));writer.writeheader();writer.writerows(samples)

def save(fig,name):
    fig.savefig(out/(name+'.png'),dpi=190)
    fig.savefig(out/(name+'.pdf'))
    plt.close(fig)

# Same physical coordinate range and color normalization for both geometries.
norm=Normalize(0,float(peaks[:,2].max()))
fig,axes=plt.subplots(1,2,figsize=(12.8,5.5),layout='constrained')
for ax,label,r in zip(axes,labels,results):
    for v in g['vv']:ax.fill(*np.array(v['pts']).T,color='#eeeeee',zorder=0)
    ax.plot(*r['wall'].T,color='#a8a8a8',lw=1.3,zorder=1)
    for n in ['floor','shelf']:ax.plot(*r['targets'][n].T,color='#555555',lw=2,zorder=2)
    h=r['maps'][2];mask=(h.segment>=0)&np.isfinite(h.q)
    ax.scatter(h.points[mask,0],h.points[mask,1],c=h.q[mask],norm=norm,cmap='inferno',s=18,zorder=3)
    idx=int(np.argmax(h.q))
    ax.plot(*h.points[idx],'o',mfc='none',mec='#2166ac',ms=12,mew=1.5,zorder=4)
    ax.set(xlim=(1.10,1.78),ylim=(-1.455,-1.19),xlabel='R (m)',ylabel='Z (m)',
           title=f'{label}  |  전체 벽 최대 {max(h.q):.3f} m^-2')
    ax.set_aspect('equal')
    ax.text(1.48,-1.228,'shelf',color='#555555',ha='center')
    ax.text(1.25,-1.408,'바닥 타깃',color='#555555',ha='center')
fig.colorbar(ScalarMappable(norm=norm,cmap='inferno'),ax=axes,label='정규화 열유속 q / P_outer (m^-2)',shrink=.70)
fig.suptitle('디버터 열유속 공간 분포 비교 · λq = 3 mm',fontsize=17)
fig.supxlabel('같은 좌표·색상 범위 / 점 = 첫 충돌 열유속 표본 / 빈 원 = 각 형상의 최대값',fontsize=10)
save(fig,'heat_spatial_3mm')

# Each column compares the same width; split shelf from floor so its smaller heat is visible.
fig,axes=plt.subplots(2,5,figsize=(17,7.3),layout='constrained')
for j,mm in enumerate(WIDTHS_MM):
    for k,n in enumerate(['floor','shelf']):
        ax=axes[k,j]
        for label,color,style,r in zip(labels,colors,['--','-'],results):
            h=r['maps'][j];mask=np.isin(h.segment,r['surface_segments'][n]);order=np.argsort(h.points[mask,0])
            ax.plot(h.points[mask,0][order],h.q[mask][order],color=color,ls=style,lw=2,label=label)
        ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.18)
        ax.set_title(f'λq = {mm} mm' if k==0 else 'shelf')
        ax.set_xlabel('첫 충돌점 R (m)')
        if j==0:ax.set_ylabel(('바닥 타깃' if k==0 else 'shelf')+'\nq / P_outer (m^-2)')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3),useMathText=False)
axes[0,0].legend(frameon=False)
fig.suptitle('표면별 열유속 분포 · 기준과 최적화 형상',fontsize=18)
fig.supxlabel('동일 패널에서 같은 축으로 비교 / 패널 간 세로축 범위는 다름 / 고정 표면은 다음 전체 벽 최대값 비교에 포함',fontsize=11)
save(fig,'heat_profiles_by_surface')

fig,axes=plt.subplots(1,2,figsize=(12,5.1),layout='constrained')
x=np.arange(5);width=.33
for i,(label,color) in enumerate(zip(labels,colors)):
    axes[0].bar(x+(i-.5)*width,peaks[i],width,color=color,label=label)
axes[0].set(xticks=x,xticklabels=WIDTHS_MM,xlabel='λq (mm)',ylabel='最大 q / P_outer (m^-2)'.replace('最大','전체 벽 최대'),title='전체 벽 최대 열유속')
axes[0].legend(frameon=False);axes[0].grid(axis='y',alpha=.18)
reduction=100*(1-peaks[1]/peaks[0])
axes[1].bar(x,reduction,color=colors[1],width=.6)
axes[1].set(xticks=x,xticklabels=WIDTHS_MM,xlabel='λq (mm)',ylabel='기준 대비 감소율 (%)',ylim=(0,50),title='최적화 형상의 최대값 감소율')
for j,v in enumerate(reduction):axes[1].text(j,v+1,f'{v:.1f}%',ha='center',fontweight='bold')
axes[1].grid(axis='y',alpha=.18)
fig.suptitle('고정 EFIT 1400 ms · 같은 outer branch 파워에서 비교',fontsize=16)
fig.supxlabel('정규화 계산 결과: 실제 MW/m² 값은 입력 파워가 있어야 환산할 수 있음',fontsize=10)
save(fig,'heat_peak_comparison')
(out/'README.md').write_text('공간 분포: heat_spatial_3mm.png\n표면별 폭 스캔: heat_profiles_by_surface.png\n전체 벽 최대값: heat_peak_comparison.png\n각 그림의 PDF 및 점별 heat_samples.csv 포함.\n공간 분포는 첫 충돌 표본이며 2D 영역 전체의 연속 열유속장이 아닙니다. q/P_outer 단위는 m^-2입니다.\n',encoding='utf-8')
print('Visualizations saved:',out,flush=True)
print('Verified peak reductions:',reduction.tolist(),flush=True)
