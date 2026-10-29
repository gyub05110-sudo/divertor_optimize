"""Render verified saved geometry without repeating optimization."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=Path(__file__).resolve().parent
a=np.load(p/'case_arrays.npz')
g=json.loads((p/'DIIID_geom.json').read_text())
rows=list(csv.DictReader((p/'width_scan.csv').open()))
peak=next(r for r in rows if r['geometry']=='candidate' and r['lambda_q_mm']=='3')
fig,axes=plt.subplots(1,2,figsize=(12,7),constrained_layout=True)
for ax in axes:
    for v in g['vv']:
        ax.fill(*np.asarray(v['pts']).T,color='.85',zorder=0)
    ax.plot(*a['baseline_wall'].T,color='.65',label='Public reference contour')
    for n,color in [('floor','tab:red'),('shelf','tab:purple')]:
        ax.plot(*a['baseline_'+n].T,'--',color=color,lw=2,label='Baseline '+n)
        ax.plot(*a['candidate_'+n].T,color=color,lw=2.5,label='Candidate '+n)
    lcfs=a['lcfs'];ax.plot(*np.vstack([lcfs,lcfs[0]]).T,color='tab:blue',lw=1.4,label='Fixed EFIT LCFS')
    ax.set(xlabel='R (m)',ylabel='Z (m)');ax.set_aspect('equal')
axes[0].set(title='Fixed equilibrium / two design surfaces')
ax=axes[1]
ax.set(xlim=(1.,1.82),ylim=(-1.49,-1.10),title='Floor: no angle limit / Shelf: 1.5 deg minimum')
ax.plot(float(peak['peak_R_m']),float(peak['peak_Z_m']),'*',color='black',ms=10,label='Candidate peak at lambda=3 mm')
ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.20),ncol=2)
fig.savefig(p/'geometry.png',dpi=180)
plt.close(fig)
