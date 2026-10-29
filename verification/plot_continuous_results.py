"""Plot saved, audited geometry and source-coordinate heat profiles."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import d3d_optimization
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--omp-folder',type=Path,help='Optional separately validated incumbent folder')
    args=parser.parse_args()
    review=json.loads((ROOT/'baseline_review/baseline_review.json').read_text(encoding='utf-8'))
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for col,source in enumerate(['IMP','OMP']):
        folder=args.omp_folder if source=='OMP' and args.omp_folder else args.output/source
        m=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        v=json.loads((folder/'viewer_result.json').read_text(encoding='utf-8'))
        data=np.load(folder/'case_arrays.npz')
        ax=axes[0,col]; base=data['baseline_wall'];candidate=data['candidate_wall']
        ax.plot(*base.T,color='#93a1ac',lw=2,label='Baseline wall')
        ax.plot(*candidate.T,color='#6550a0',lw=1.5,label='Selected candidate')
        for name,color in [('inner','#bc35a2'),('outer','#d34730')]:
            leg=np.asarray(review['legs'][name]['points'])
            ax.plot(*leg.T,color=color,lw=1.5,label=name+' separatrix')
        for surface in m['surfaces']:
            if surface['original_segment'] in [69,70]:
                points=np.asarray(surface['endpoints_RZ_m'])
                ax.plot(*points.T,color='#008d99',lw=3,alpha=.6)
        w=v.get('boundary_witness')
        if w and w['resolved']:ax.scatter(*w['point'],color='black',s=40,marker='x',zorder=5,label='SOL limit witness')
        ax.axhline(-1.15,color='#008d99',ls=':',lw=1)
        ax.set(xlim=(.99,1.72),ylim=(-1.41,-1.1),xlabel='R [m]',ylabel='Z [m]',aspect='equal')
        status=m['validity_audit_summary'];ratio=status.get('finest_worst_peak_ratio',float('nan'))
        origin='retained incumbent' if m.get('replay_only') else 'new search'
        ax.set_title(f'{source} ({origin}): ratio {ratio:.6f}; reduction {100*(1-ratio):.3f}%')
        ax.legend(fontsize=7,loc='lower right');ax.grid(alpha=.2)
        ax=axes[1,col]
        qbase=data['baseline_q_3mm'];qcand=data['candidate_q_3mm']
        edges=m['support_mm']*np.linspace(0,1,len(qbase)+1)**2
        x=(edges[:-1]+edges[1:])/2
        bq=m['baseline']['boundary_witness']['q'][2]
        cq=w['q'][2]
        ax.plot(np.r_[0,x],np.r_[bq,qbase],color='#8796a2',label='Baseline')
        ax.plot(np.r_[0,x],np.r_[cq,qcand],color='#6550a0',label='Candidate')
        ax.scatter([0,0],[bq,cq],color=['#8796a2','#6550a0'],marker='x',s=35,zorder=4)
        ax.set(xlim=(0,20),ylim=(0,None),xlabel=f'{source} upstream radial offset [mm]',ylabel='First-hit q/P [1/m²]',
               title='λq = 3 mm; whole-wall first hit + x→0 limit')
        ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Continuous-boundary numerical audit | fixed equilibrium, no cross-field transport\nIndependent unit-power sources; no experimental validation or global-optimality certificate',fontsize=12)
    fig.savefig(args.output/'boundary_results.png',dpi=180)
    fig.savefig(args.output/'boundary_results.pdf')
    plt.close(fig)

if __name__=='__main__':main()
