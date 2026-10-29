"""Read-only, independently reconstructed routes excluded by the legacy kernel."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import json,csv
import numpy as np
import shapely
from shapely.geometry import Polygon
from d3d_optimization.broadening import _cdf_interval


def legacy_routes(paths,hits,dep,weights,lcfs,solids,anchor_m):
    w=np.asarray(weights);mesh=dep['wall'];xi=dep['cell_xi_m']
    domain=Polygon(mesh).buffer(1e-10)
    core=None if lcfs is None else (lcfs if hasattr(lcfs,'geom_type') else Polygon(lcfs))
    solid=shapely.union_all([solids] if hasattr(solids,'geom_type') else solids)
    rows=[]
    for i,weight in enumerate(w):
        r=dict(donor=i,weight=float(weight),branch_tail=0.,blocked_lcfs=0.,blocked_solid=0.,blocked_domain=0.,deposited=0.,unhit=0.)
        if hits.segment[i]<0:r['unhit']=float(weight);rows.append(r);continue
        sigma=dep['sigma_m'][i];center=dep['source_xi_m'][i]
        branch=dep['diagnostics']['donor_branch'][i]
        cells=np.flatnonzero(dep['cell_branch']==branch)
        lo=(xi[cells].min(axis=1)-center)/sigma;hi=(xi[cells].max(axis=1)-center)/sigma
        cells=cells[(hi>-8)&(lo<8)&(hi>lo)]
        mass=_cdf_interval((xi[cells].min(axis=1)-center)/sigma,(xi[cells].max(axis=1)-center)/sigma)
        r['branch_tail']=float(weight*max(0,1-mass.sum()))
        anchor=np.asarray(paths[i].interpolate(hits.distance[i]-anchor_m).coords[0])
        nodes=np.stack([np.broadcast_to(anchor,mesh[cells].shape),mesh[cells],mesh[cells+1]],axis=1)
        fan=shapely.polygons(nodes);zero=shapely.area(fan)==0
        if np.any(zero):fan[zero]=shapely.linestrings(nodes[zero])
        core_block=np.zeros(len(cells),bool) if core is None else shapely.intersects(core,fan)
        solid_block=shapely.intersects(solid,fan)&~core_block
        domain_block=~shapely.covers(domain,fan)&~core_block&~solid_block
        visible=~(core_block|solid_block|domain_block)
        for key,mask in [('blocked_lcfs',core_block),('blocked_solid',solid_block),('blocked_domain',domain_block),('deposited',visible)]:
            r[key]=float(weight*mass[mask].sum())
        r.update(hit_R=float(hits.points[i,0]),hit_Z=float(hits.points[i,1]),sigma_m=float(sigma),source_xi_m=float(center),
                 branch_min_xi=float(xi[dep['cell_branch']==branch].min()),branch_max_xi=float(xi[dep['cell_branch']==branch].max()))
        rows.append(r)
    keys=['branch_tail','blocked_lcfs','blocked_solid','blocked_domain','deposited','unhit']
    totals={k:sum(r[k] for r in rows) for k in keys}
    totals['blocked']=sum(totals[k] for k in keys if k.startswith('blocked_'))
    totals['unresolved']=totals['blocked']+totals['branch_tail']+totals['unhit']
    error=max(abs(totals[k]-float(np.asarray(dep['ledger'][k]).ravel()[0])) for k in ['deposited','blocked','branch_tail','unhit','unresolved'])
    return dict(totals=totals,donors=rows,reconciliation_max_error=error,
                blocked_partition='exclusive: LCFS, then solid, then outside vessel; these are legacy exclusions, not physical destinations')


def main():
    import argparse
    from d3d_optimization.independent_sources import prepare_evaluator
    from d3d_optimization.reference_transport import ReferencePoint
    from d3d_optimization.transport_evaluator import TransportEvaluator
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    d,meta=prepare_evaluator(ROOT,'OMP',count=320,step=.0005,nodes=161,boundary_mode='continuous')
    e=TransportEvaluator(d,meta,ReferencePoint(),wall_step_m=.001)
    result=e.baseline;dep=result['deposition'];hits=e._hits(d.baseline)
    report=legacy_routes(e.paths,hits,dep,e.sources[0].weights,e.lcfs,e.vessel_solids,e.controls.anchor_m)
    for i,row in enumerate(report['donors']):
        row.update(x_mm=float(e.sources[0].x[i]*1000),Tt_eV=float(result['thermal']['Tt_eV'][i]),L_m=float(result['thermal']['length_m'][i]))
    report['cold_tail_unresolved']=sum(r['branch_tail']+sum(r[k] for k in ['blocked_lcfs','blocked_solid','blocked_domain']) for r in report['donors'] if r['Tt_eV']<5)
    report['metadata']=meta
    (args.output/'diagnostic.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf8')
    with (args.output/'donors.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(report['donors'][0]));writer.writeheader();writer.writerows(report['donors'])
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    lost=np.array([r['branch_tail']+r['blocked_lcfs']+r['blocked_solid']+r['blocked_domain'] for r in report['donors']])
    ax[0].plot(d.baseline['wall'][:,0],d.baseline['wall'][:,1],color='gray')
    keep=lost>1e-7
    s=ax[0].scatter(hits.points[keep,0],hits.points[keep,1],c=lost[keep]*100,s=15+5000*lost[keep],cmap='magma')
    fig.colorbar(s,ax=ax[0],label='Excluded power per donor (% input)')
    ax[0].set(xlim=(1.1,1.6),ylim=(-1.42,-1.2),xlabel='R (m)',ylabel='Z (m)',title='Legacy exclusions at fixed first hits')
    for key in ['branch_tail','blocked_lcfs','blocked_solid','blocked_domain']:
        ax[1].plot(e.sources[0].x*1000,[r[key]*100 for r in report['donors']],label=key)
    ax[1].set(xlim=(0,25),xlabel='Source distance (mm)',ylabel='Excluded input power per donor (%)',title='Cause of excluded power')
    ax[1].legend();fig.tight_layout();fig.savefig(args.output/'diagnostic.png',dpi=160);plt.close(fig)
    print(json.dumps({k:v for k,v in report.items() if k not in ['donors','metadata']},indent=2))


if __name__=='__main__':main()
