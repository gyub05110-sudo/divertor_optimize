"""Inspect saved multi-peak results; does not change any qualification criterion."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qualify_omp_stages import read,checked_read,save,sha
from export_qualified_finals import STUDY


def main():
    identity='0eb363d2b21c1d391948f380244813840964a06edab792fffca167be74755cc1'
    folder=STUDY/'R020/validation'/identity;v=checked_read(folder/'validation.json')
    if not v['complete']:raise ValueError('wait for complete validation')
    out=STUDY/'diagnostics';out.mkdir(exist_ok=True)
    fine=v['levels'][2]['rows']['selected'];quintic=v['levels'][3]['rows']['selected'];points=[fine['peak_RZ_m'],quintic['peak_RZ_m']]
    rows=[];fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for level,color in zip(v['levels'],['#aaa','#5686aa','#135a78','#b83f6c']):
        name=f"{level['method']}_{level['tubes']}";path=folder/(name+'.npz')
        if sha(path)!=v['artifacts'][path.name]:raise ValueError('saved arrays changed')
        with np.load(path) as z:
            p=z['points'];q=z['q_MW_m2'];mesh=z['receiving_wall'];centers=.5*(mesh[:-1]+mesh[1:])
            mean_q=z['cell_power']/np.diff(z['area_edges']) # assigned power = 1 MW
            mask=(centers[:,0]>1.17)&(centers[:,0]<1.43)&(centers[:,1]<-1.34)
            ax.plot(centers[mask,0],mean_q[mask],color=color,label=name.replace('_',' '),lw=1.6)
            near=[];local=[]
            for target in points:
                index=int(np.argmin(np.linalg.norm(p-target,axis=1)))
                near.append(dict(target_RZ=target,sample_RZ=p[index].tolist(),q_MW_m2=float(q[index]),distance_m=float(np.linalg.norm(p[index]-target))))
                ids=np.flatnonzero(np.linalg.norm(p-target,axis=1)<=.001)
                local_index=int(ids[np.argmax(q[ids])])
                local.append(dict(neighborhood_radius_m=.001,sample_RZ=p[local_index].tolist(),q_MW_m2=float(q[local_index]),distance_m=float(np.linalg.norm(p[local_index]-target))))
        rows.append(dict(level=name,maximum_RZ=level['rows']['selected']['peak_RZ_m'],maximum_q=level['rows']['selected']['peak_MW_m2'],near_competing_peaks=near,stored_neighborhood_maxima=local))
    for i,p in enumerate(points):ax.axvline(p[0],color=['#135a78','#b83f6c'][i],ls='--',alpha=.6)
    ax.set(xlim=(1.20,1.31),ylim=(1.42,1.45),xlabel='R along the lower OMP target [m]',ylabel='Cell-average q [MW/m²]',
        title='OMP 20 mm, archived stage 18: near-peak zoom\nArgmax changes branch; the full footprint does not translate by 71 mm')
    ax.grid(alpha=.2);ax.legend();fig.savefig(out/'omp20-competing-peaks.png',dpi=180);plt.close(fig)
    save(out/'peak-location.json',dict(candidate=identity,criteria_unchanged=True,verdict=v['verdict'],levels=rows,
        generator_sha256=sha(Path(__file__)),validation_sha256=sha(folder/'validation.json'),
        interpretation='Both near-equal local peaks persist. Interpolation changes their ordering. This is an argmax-location failure, not evidence of a rigid footprint translation or failed geometry. Small CDF/value changes do not establish a stable hottest-point location or uniform pointwise convergence. Existing 2 mm criterion remains enforced.'))
    (out/'diagnose_peak_location.py').write_bytes(Path(__file__).read_bytes())
    print(out/'peak-location.json')


if __name__=='__main__':main()
