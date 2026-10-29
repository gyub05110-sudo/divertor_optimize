"""Render the inspected equilibrium and sampled first-hit regions for review."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'baseline_review'
d=json.loads((OUT/'baseline_review.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False,'font.size':10,'savefig.dpi':180})
colors={'IMP':'#008d99','OMP':'#36a25b','inner':'#bc35a2','outer':'#d34730','design':'#ed9c18'}
wall=np.array(d['wall']);lcfs=np.array(d['lcfs']);xp=np.array(d['xpoint_RZ_m']);pfr=np.array(d['pfr_polygon'])
fig,axs=plt.subplots(1,2,figsize=(15,8),gridspec_kw={'width_ratios':[1,1.15]})

def draw(ax,detail=False):
    ax.plot(*wall.T,color='#738392',lw=1.8,zorder=3)
    ax.plot(*lcfs.T,color='#4e77ad',lw=1.1,alpha=.65,zorder=2)
    ax.fill(*pfr.T,facecolor='#efcfe1',edgecolor='#b875a2',alpha=.65,hatch='///',zorder=1)
    for seg in d['existing_design_segments']:
        ax.plot(*wall[seg:seg+2].T,color=colors['design'],lw=7,alpha=.9,solid_capstyle='round',zorder=4)
    for row in d['footprints']:
        p=np.array(row['points']);ax.plot(*p.T,color='#9cabb5',lw=2.4,zorder=4,solid_capstyle='round')
    for row in d['optimization_footprints']:
        p=np.array(row['points']);ax.plot(*p.T,color=colors[row['source']],lw=4,zorder=5,solid_capstyle='round')
    for source,s in d['sources'].items():
        for row in s['paths']:
            p=np.array(row['points']);ax.plot(*p.T,color=colors[source],lw=.8,alpha=.5,zorder=2)
    for name,leg in d['legs'].items():
        p=np.array(leg['points']);ax.plot(*p.T,color=colors[name],lw=2.5,zorder=6)
        end=np.array(leg['strike_RZ_m']);ax.scatter(*end,s=30,color=colors[name],zorder=8)
    ax.scatter(*xp,s=80,marker='x',color='#131e30',lw=2.2,zorder=8)
    ax.set_aspect('equal');ax.set_xlabel('R (m)');ax.set_ylabel('Z (m)')
    ax.grid(alpha=.13)
    if detail:
        ax.set_xlim(.99,1.78);ax.set_ylim(-1.405,-.82)
        ax.annotate('X-point\n(1.17659, -1.25240)',xy=xp,xytext=(1.27,-1.09),arrowprops={'arrowstyle':'-','color':'#24354b'},fontsize=10)
        ax.annotate('inner leg',xy=d['legs']['inner']['points'][int(len(d['legs']['inner']['points'])*.6)],xytext=(1.04,-1.11),color=colors['inner'],arrowprops={'arrowstyle':'-','color':colors['inner']})
        ax.annotate('outer leg',xy=d['legs']['outer']['points'][int(len(d['legs']['outer']['points'])*.6)],xytext=(1.30,-1.29),color=colors['outer'],arrowprops={'arrowstyle':'-','color':colors['outer']})
        ax.annotate('PFR\n이 모델로 열분포 미정',xy=(1.161,-1.327),xytext=(1.41,-1.397),arrowprops={'arrowstyle':'-','color':'#95547f'},color='#854269',fontsize=10)
        ax.text(1.36,-1.335,'floor · seg 68',color='#966116',ha='left')
        ax.text(1.50,-1.222,'shelf · seg 59',color='#966116')
        ax.annotate('IMP 설계영역\nZ=-1.15 m부터 아래쪽\n+ 경사면 일부 seg 69',xy=(1.01749,-1.18),xytext=(1.28,-.90),arrowprops={'arrowstyle':'-','color':colors['IMP']},color=colors['IMP'])
        ax.plot([1.005,1.09],[-1.15,-1.15],ls='--',color=colors['IMP'],lw=1.2,zorder=7)
        ax.text(1.04,-1.143,'Z=-1.15',color=colors['IMP'],fontsize=9)
    else:
        ax.set_xlim(.90,2.46);ax.set_ylim(-1.45,.14)
        z=d['axis_RZ_m'][1];ax.axhline(z,color='#24354b',ls='--',lw=.8)
        for name,s in d['sources'].items():
            r,z=s['midplane_RZ_m'];ax.scatter(r,z,color=colors[name],s=45,zorder=8)
            ax.annotate(name+' 열원 위치',xy=(r,z),xytext=(r+(.11 if name=='IMP' else -.26),z+.09),color=colors[name],arrowprops={'arrowstyle':'-','color':colors[name]})
        ax.text(1.70,-.20,'중간면 Z = -0.06456 m',ha='center',color='#536577',fontsize=9)
        ax.annotate('수직 벽 상부는 고정\n설계영역은 Z<=-1.15 m',xy=(1.0174,-.52),xytext=(1.22,-.49),color='#687d8d',arrowprops={'arrowstyle':'-','color':'#687d8d'})
        ax.annotate('첨부한 경사 구간: 최적화 제외\n고정 벽으로 열부하 계산 유지',xy=(1.98,-.983),xytext=(1.50,-.73),color='#687d8d',arrowprops={'arrowstyle':'-','color':'#687d8d'})
        ax.add_patch(plt.Rectangle((.99,-1.405),.79,.585,fill=False,ls=':',ec='#8a97a6'))

draw(axs[0]);draw(axs[1],True)
axs[0].set_title('① 실제 IMP·OMP 위치와 최초 충돌 경로',loc='left',fontweight='bold',pad=15)
axs[1].set_title('② 두 separatrix leg와 디버터 설계영역',loc='left',fontweight='bold',pad=15)
fig.suptitle('사진의 기준 형상에서 새 요청사항 확인 · DIII-D / EFIT 1400 ms',fontsize=17,fontweight='bold',x=.08,ha='left',y=.98)
legend=[Line2D([0],[0],color=colors['design'],lw=6,label='기존 설계영역: floor + shelf 전체'),
        Line2D([0],[0],color=colors['IMP'],lw=4,label='선택 IMP 설계영역: Z<=-1.15 m'),
        Line2D([0],[0],color=colors['OMP'],lw=4,label='선택 OMP 설계영역: 경사 벽 일부 제외'),
        Line2D([0],[0],color=colors['inner'],lw=2.5,label='inner separatrix leg'),
        Line2D([0],[0],color=colors['outer'],lw=2.5,label='outer separatrix leg'),
        Patch(facecolor='#efcfe1',hatch='///',label='PFR: 직접 SOL 매핑으로 열분포 미정')]
fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.055),ncol=3,frameon=False,fontsize=10)
fig.text(.08,.02,'색은 열유속이 아니라 경로·영역 구분입니다. 후보 최적화는 아직 수행하지 않았으며, 충돌 표본 범위의 연속 경계는 이후 세분화합니다.',fontsize=10,color='#526577')
fig.subplots_adjust(left=.07,right=.98,top=.89,bottom=.18,wspace=.20)
fig.savefig(OUT/'baseline_requested_regions.png')
fig.savefig(OUT/'baseline_requested_regions.pdf')
plt.close(fig)
print(OUT/'baseline_requested_regions.png')
