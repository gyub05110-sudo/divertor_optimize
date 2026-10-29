"""Render sampled heat as colored lines on the complete saved wall geometry."""
import csv
import json
from pathlib import Path
import numpy as np
from shapely.geometry import LineString
from shapely.ops import substring
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle

p=Path(__file__).resolve().parent
out=p/'comparison'
a=np.load(p/'case_arrays.npz')
g=json.loads((p/'DIIID_geom.json').read_text())
data=list(csv.DictReader((out/'heat_samples.csv').open(encoding='utf-8-sig')))
summary=list(csv.DictReader((p/'width_scan.csv').open()))
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False,'font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False})
labels=['기준 형상','최적화 형상']
keys=['baseline','candidate']

def colored_wall(wall,rows):
    line=LineString(wall)
    samples=[]
    for row in rows:
        if row['component']=='미포착':continue
        point=np.array([float(row['R_m']),float(row['Z_m'])])
        from shapely.geometry import Point
        assert line.distance(Point(point))<1e-7
        samples.append((line.project(Point(point)),float(row['q_over_P_outer_per_m2']),row['component']))
    samples.sort()
    pieces=[];values=[]
    for (s0,q0,c0),(s1,q1,c1) in zip(samples[:-1],samples[1:]):
        # Join only local samples on the same component. Unsampled wall stays grey.
        if s1<=s0 or s1-s0>.02 or c0!=c1:continue
        ss=np.linspace(s0,s1,9);qq=np.linspace(q0,q1,9)
        for j in range(8):
            piece=substring(line,ss[j],ss[j+1])
            pieces.append(np.asarray(piece.coords));values.append((qq[j]+qq[j+1])/2)
    return pieces,np.array(values)

for mm in [3,10]:
    selections=[[r for r in data if r['geometry']==label and int(r['lambda_q_mm'])==mm] for label in labels]
    peaks=[max(float(r['q_over_P_outer_per_m2']) for r in rows) for rows in selections]
    for key,peak in zip(keys,peaks):
        recorded=next(float(r['all_wall_peak_q_over_P']) for r in summary if r['geometry']==key and int(r['lambda_q_mm'])==mm)
        assert abs(recorded-peak)<1e-10
    norm=Normalize(0,max(peaks));cmap='inferno'
    fig,axes=plt.subplots(2,2,figsize=(12.5,11.2),layout='constrained',gridspec_kw={'height_ratios':[1.7,1]})
    zoom=(1.10,1.78,-1.46,-1.18)
    for j,(key,label,rows,peak) in enumerate(zip(keys,labels,selections,peaks)):
        wall=a[key+'_wall'];pieces,values=colored_wall(wall,rows)
        peakrow=max(rows,key=lambda r:float(r['q_over_P_outer_per_m2']))
        peakpoint=[float(peakrow['R_m']),float(peakrow['Z_m'])]
        for k in range(2):
            ax=axes[k,j]
            for v in g['vv']:ax.fill(*np.array(v['pts']).T,color='#eeeeee',zorder=0)
            ax.plot(*wall.T,color='#aaaaaa',lw=1.5,zorder=1)
            lcfs=a['lcfs'];ax.plot(*np.vstack([lcfs,lcfs[0]]).T,color='#3976a8',lw=1.1,ls='--',zorder=1)
            ax.add_collection(LineCollection(pieces,colors='#454545',linewidths=8,capstyle='round',zorder=2))
            colored=LineCollection(pieces,cmap=cmap,norm=norm,linewidths=6,capstyle='round',zorder=3)
            colored.set_array(values);ax.add_collection(colored)
            ax.plot(*peakpoint,'o',mfc='none',mec='#00a6b2',ms=13,mew=2,zorder=5)
            ax.set_aspect('equal');ax.set(xlabel='R (m)',ylabel='Z (m)')
            if k==0:
                ax.set(xlim=(.85,2.55),ylim=(-1.58,1.58),title=f'{label} · 최대 {peak:.3f} m^-2')
                ax.add_patch(Rectangle((zoom[0],zoom[2]),zoom[1]-zoom[0],zoom[3]-zoom[2],
                                      fill=False,ec='#666666',lw=1,ls=':',zorder=4))
                ax.text(1.7,.4,'고정 평형 경계',color='#3976a8',fontsize=10)
            else:
                ax.set(xlim=zoom[:2],ylim=zoom[2:],title='하부 디버터 확대')
                ax.text(1.55,-1.218,'shelf',ha='center',color='#555555')
                ax.text(1.24,-1.420,'바닥 타깃',ha='center',color='#555555')
                ax.annotate(f'최대 {peak:.3f}',xy=peakpoint,xytext=(1.12,-1.205),fontsize=10,
                            arrowprops={'arrowstyle':'->','color':'#00a6b2'},color='#00808a')
    fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),shrink=.75,
                 label='정규화 열유속 q / P_outer (m^-2) · 두 형상 공통 범위',pad=.025)
    reduction=100*(1-peaks[1]/peaks[0])
    fig.suptitle(f'전체 DIII-D 2D 형상의 표면 열유속 · λq = {mm} mm\n최적화 형상의 전체 벽 최대값 {reduction:.1f}% 감소',fontsize=17)
    fig.supxlabel('색 선 = 첫 충돌 표본 사이의 표면 열유속 표시 / 회색 = 표본이 없는 벽 / 청록색 원 = 최대값\n선 두께는 가독성을 위해 일정하게 표시 · outer branch 직접 열유속 모델 · 실제 MW/m²가 아님',fontsize=10)
    name=f'full_geometry_heat_{mm}mm'
    fig.savefig(out/(name+'.png'),dpi=190)
    fig.savefig(out/(name+'.pdf'))
    plt.close(fig)
    print('Saved',name,'peak reduction',reduction,flush=True)
