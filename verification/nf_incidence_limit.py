"""Adversarial limit: radial displacement before first interception of a smooth wall.
The original connected-wall model and the analytic flux-coordinate reference differ in
where transverse transport acts. A discrepancy diagnoses model form in this limit.
"""
import argparse
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from d3d_optimization.flux_footprint import branch_power
from d3d_optimization.boundary_transport import deposit_boundary
from d3d_optimization.frozen_mapping import Hits
from d3d_optimization.broadening_study import save


class RadialLabel:
    @staticmethod
    def gradient(points):
        return np.tile([1.,0.],(len(points),1))


def compare_target(slope,curvature,sigma):
    # Flux label x = R-R0. Field lines are vertical in the poloidal section.
    # A single upstream Gaussian packet has an exactly known CDF on the target.
    r0=2.
    x=np.linspace(-.04,.04,641)
    target=np.column_stack((r0+x,slope*x+.5*curvature*x*x))
    top=float(target[:,1].max()+.06)
    wall=np.vstack([target,[target[-1,0],top],[target[0,0],top],target[0]])
    hits=Hits(np.array([len(x)//2]),np.array([[r0,0.]]),np.array([1.]))
    d=deposit_boundary(RadialLabel(),hits,wall,[1.],[sigma],1.,.00025)
    cells=d['cell_parent']<len(x)-1
    radial_edges=d['wall'][:cells.sum()+1,0]-r0
    ref=branch_power(radial_edges,[0.],[sigma],[1.])
    physical=np.zeros_like(d['cell_power'])
    physical[cells]=ref['cell_power']
    q=physical/d['areas']
    orig_q=d['cell_power']/d['areas']
    # Compare against analytic density at segment centers, separately from the CDF integrator.
    mid=.5*(radial_edges[1:]+radial_edges[:-1])
    grad=slope+curvature*mid
    exact=np.exp(-.5*(mid/sigma)**2)/(np.sqrt(2*np.pi)*sigma*2*np.pi*(r0+mid)*np.sqrt(1+grad*grad))
    return dict(slope=slope,curvature_per_m=curvature,sigma_m=sigma,
                reference_peak=float(q.max()),wall_kernel_peak=float(orig_q.max()),
                analytic_midpoint_peak=float(exact.max()),ref_peak_cell_relative=abs(float(q.max()/exact.max()-1)),
                reference_A2=float(physical.sum()**2/np.sum(physical**2/d['areas'])),
                wall_kernel_A2=float(d['cell_power'].sum()**2/np.sum(d['cell_power']**2/d['areas'])),
                total_variation=.5*float(np.abs(physical-d['cell_power']).sum()+ref['below']+ref['above']),
                reference_unresolved=ref['below']+ref['above'],original_balance=d['ledger']['balance_error'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    rows=[compare_target(slope,kappa,.002) for slope,kappa in ((0.,0.),(.5,0.),(1.,0.),(2.,0.),(4.,0.),(2.,120.))]
    save(a.output,{'source_power':1.,'rows':rows,'scope':'single donor, monotone radial label, pre-impact spreading'},exclusive=True)
    for row in rows:
        print(row['slope'],row['curvature_per_m'],row['reference_peak'],row['wall_kernel_peak'],row['total_variation'])

if __name__=='__main__':main()

