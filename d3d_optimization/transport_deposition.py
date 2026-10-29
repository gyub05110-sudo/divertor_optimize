"""Explicit donor-width version of the reviewed endpoint kernel.

Copied from broadening.py on 2026-10-28 to preserve historical calculations.
The only numerical change is accepting sigma_m instead of deriving sigma(L).
All branch, visibility, CDF and area-density algorithms remain identical.
"""
from dataclasses import dataclass, asdict
import numpy as np
import shapely
from shapely.geometry import LineString, Polygon
from scipy.optimize import brentq
from .broadening import _wall_mesh, _cdf_interval


@dataclass(frozen=True)
class DepositionControls:
    wall_step_m: float = .002
    anchor_m: float = .005
    missing_tolerance: float = .001

    def __post_init__(self):
        if (not all(np.isfinite(v) for v in asdict(self).values()) or
                self.wall_step_m<=0 or self.anchor_m<=0 or not 0<=self.missing_tolerance<1):
            raise ValueError('invalid deposition controls')


def deposit_sigma(field, paths, lengths, hits, wall, weights_matrix,
                     source_rz, source_sign, lcfs, vessel_solids, config, sigma_m):
    """Deposit all source widths in a single geometry/kernel pass.

    `source_rz` is the separatrix point. Actual donor flux is sampled at each
    path start. Weights have shape (donor, source width). Parent segment IDs
    refer to supplied candidate wall segments, not original design segments.
    """
    w=np.asarray(weights_matrix,float)
    if w.ndim==1: w=w[:,None]
    n=len(paths)
    if w.shape[0]!=n or not np.all(np.isfinite(w)) or np.any(w<0) or np.any(w.sum(axis=0)<=0):
        raise ValueError('finite nonnegative N by W source powers required')
    sep=np.asarray(source_rz,float).reshape(1,2)
    grad=float(field.gradient(sep)[0,0]); offset=float(field.values(sep)[0])
    if source_sign not in (-1,1) or not np.isfinite(grad) or grad==0:
        raise ValueError('source sign and nonzero separatrix radial flux gradient required')
    scale=np.sign(grad*source_sign)/abs(grad)
    mesh,parent,ds,tangent,xi,branch=_wall_mesh(field,wall,scale,offset,config.wall_step_m)
    m=len(parent); wall=np.asarray(wall,float)
    vessel=Polygon(wall)
    if not vessel.is_valid or vessel.area<=0: raise ValueError('invalid vessel polygon')
    # Tolerance only compensates GEOS roundoff on material endpoints.
    domain=vessel.buffer(1e-10)
    obstacles=[]
    solids=([vessel_solids] if hasattr(vessel_solids,'geom_type') else
            ([] if vessel_solids is None else list(vessel_solids)))
    for geom in ([lcfs] if lcfs is not None else [])+solids:
        poly=geom if hasattr(geom,'geom_type') else Polygon(geom)
        if not poly.is_valid: raise ValueError('invalid material/LCFS geometry')
        obstacles.append(poly)
    obstacle=shapely.union_all(obstacles) if obstacles else None
    shapely.prepare(domain)
    if obstacle is not None: shapely.prepare(obstacle)
    active=w.sum(axis=1)>0
    hit=np.asarray(hits.segment)>=0
    if np.asarray(hits.segment).shape!=(n,) or np.any(np.asarray(hits.segment)[hit]>=len(wall)-1):
        raise ValueError('invalid first-hit segment indices')
    distance=np.where(active&hit,np.asarray(hits.distance),np.inf)
    L=lengths.at(distance); sigma=np.asarray(sigma_m,float).copy()
    if sigma.shape!=(n,): raise ValueError("one explicit sigma per donor required")
    good=active&hit
    if np.any(~np.isfinite(sigma[good])) or np.any(sigma[good]<=0):
        raise ValueError('positive finite donor connection lengths required')
    source_xi=scale*(np.array([p[0] for p in lengths.psi])-offset)
    donor_cell=np.full(n,-1,int); donor_branch=np.full(n,-1,int)
    anchors=np.full((n,2),np.nan); drift=np.full(n,np.nan)
    for i in np.flatnonzero(good):
        candidates=np.flatnonzero(parent==hits.segment[i])
        point=np.asarray(hits.points[i]); a=mesh[candidates]; d=mesh[candidates+1]-a
        f=np.clip(((point-a)*d).sum(axis=1)/(d*d).sum(axis=1),0,1)
        residual=np.linalg.norm(a+f[:,None]*d-point,axis=1)
        j=candidates[np.argmin(residual)]
        if np.min(residual)>1e-7: raise ValueError('first hit not on indicated wall segment')
        # At a physical branch endpoint both adjacent branches are possible;
        # reject this ambiguity rather than select whichever has the lower peak.
        vertices=np.flatnonzero(np.linalg.norm(mesh-point,axis=1)<1e-10)
        for vertex in vertices:
            left=(vertex-1)%m; right=vertex%m
            if branch[left]!=branch[right]:
                raise ValueError('ambiguous first hit at receiving branch boundary')
        donor_cell[i]=j; donor_branch[i]=branch[j]
        path=paths[i] if hasattr(paths[i],'interpolate') else LineString(paths[i])
        if distance[i]<config.anchor_m: raise ValueError('first hit too close for fixed upstream anchor')
        anchors[i]=np.asarray(path.interpolate(distance[i]-config.anchor_m).coords[0])
        before=lengths.s[i]<=distance[i]
        drift[i]=abs(scale)*max(float(np.max(np.abs(lengths.psi[i][before]-lengths.psi[i][0]))),
                                  abs(float(field.values([point])[0])-lengths.psi[i][0]))
    # Five points per cell and donor-centered flux roots capture narrow peaks.
    fractions=np.array([0.,.25,.5,.75,1.])
    points=(mesh[:-1,None,:]+fractions[None,:,None]*(mesh[1:]-mesh[:-1])[:,None,:]).reshape(-1,2)
    sample_cell=np.repeat(np.arange(m),len(fractions))
    extra_points=[]; extra_cells=[]
    for i in np.flatnonzero(good):
        cells=np.flatnonzero(branch==donor_branch[i])
        lo=np.minimum(xi[cells],xi[cells+1]); hi=np.maximum(xi[cells],xi[cells+1])
        for target in source_xi[i]+sigma[i]*np.array([-2.,-1.,-.5,0.,.5,1.,2.]):
            found=cells[(lo<target)&(hi>target)]
            if len(found):
                j=int(found[0]); a,b=mesh[j:j+2]
                f=brentq(lambda u: float(scale*(field.values([a+u*(b-a)])[0]-offset)-target),0,1,xtol=1e-12)
                extra_points.append(a+f*(b-a)); extra_cells.append(j)
    if extra_points:
        points=np.vstack([points,extra_points]); sample_cell=np.r_[sample_cell,extra_cells]
    sample_xi=scale*(np.asarray(field.values(points))-offset)
    sample_jac=abs((scale*np.asarray(field.gradient(points))*tangent[sample_cell]).sum(axis=1))/(2*np.pi*points[:,0])
    if not np.all(np.isfinite(sample_jac)): raise ValueError('nonfinite surface Jacobian')
    # Group point indices once; kernel support is sparse even for a dense wall.
    order=np.argsort(sample_cell,kind='stable'); counts=np.bincount(sample_cell,minlength=m)
    sample_indices=np.split(order,np.cumsum(counts)[:-1])
    q=np.zeros((len(points),w.shape[1])); power=np.zeros((m,w.shape[1]))
    deposited=np.zeros(n); blocked=np.zeros(n); tail=np.zeros(n)
    visibility_tests=0
    for i in np.flatnonzero(good):
        lo=(np.minimum(xi[:-1],xi[1:])-source_xi[i])/sigma[i]
        hi=(np.maximum(xi[:-1],xi[1:])-source_xi[i])/sigma[i]
        cells=np.flatnonzero((branch==donor_branch[i])&(hi>-8)&(lo<8)&(hi>lo))
        mass=_cdf_interval(lo[cells],hi[cells])
        # Full selected-cell mass is evaluated (including its tiny far tail).
        # Every excluded cell/tail is explicitly included in branch_tail.
        total=float(mass.sum())
        if total>1+1e-11 or np.any(mass<0): raise ValueError('invalid Gaussian interval mass')
        tail[i]=max(0.,1-total)
        if not len(cells): continue
        # The full triangular fan is the union of ALL anchor-to-cell chords.
        # Point probes can miss a thin solid between their rays; checking this
        # region instead guarantees a common conservative mask for CDF and q.
        nodes=np.stack([np.broadcast_to(anchors[i],mesh[cells].shape),
                        mesh[cells],mesh[cells+1]],axis=1)
        fans=shapely.polygons(nodes)
        degenerate=shapely.area(fans)==0
        if np.any(degenerate):
            # Collinear anchor/cell: all chords form this line, not an invalid
            # zero-area polygon. Backtracking vertices retain the full span.
            fans[degenerate]=shapely.linestrings(nodes[degenerate])
        visible=np.asarray(shapely.covers(domain,fans))
        if obstacle is not None: visible &= ~np.asarray(shapely.intersects(obstacle,fans))
        visibility_tests+=len(fans)
        deposited[i]=float(mass[visible].sum()); blocked[i]=float(mass[~visible].sum())
        cells=cells[visible]; mass=mass[visible]
        power[cells]+=mass[:,None]*w[i]
        if len(cells):
            indices=np.concatenate([sample_indices[c] for c in cells])
            z=(sample_xi[indices]-source_xi[i])/sigma[i]
            density=np.exp(-.5*z*z)/(np.sqrt(2*np.pi)*sigma[i])*sample_jac[indices]
            q[indices]+=density[:,None]*w[i]
    prescribed=w.sum(axis=0)
    ledger=dict(prescribed=prescribed,deposited=deposited@w,blocked=blocked@w,
                branch_tail=tail@w,unhit=w[~hit].sum(axis=0))
    ledger['unresolved']=ledger['blocked']+ledger['branch_tail']+ledger['unhit']
    ledger['balance_error']=abs(ledger['deposited']+ledger['unresolved']-prescribed)/prescribed
    reasons=[]
    if np.any(ledger['unresolved']/prescribed>config.missing_tolerance): reasons.append('unresolved_power_exceeds_tolerance')
    if np.any(ledger['balance_error']>1e-10): reasons.append('power_balance_error')
    if not np.all(np.isfinite(q)) or np.any(q<0): raise ValueError('invalid surface heat density')
    # Independent composite Boole integration of sampled q*2*pi*R ds.
    # This deliberately retains quadrature error for the refinement audit.
    boole=np.array([7.,32.,12.,32.,7.])/90
    quadrature_area=(2*np.pi*points[:5*m,0].reshape(m,5)*ds[:,None]*boole).ravel()
    quadrature_power=quadrature_area@q[:5*m]
    return dict(quadrature_integrated_power=quadrature_power,cell_xi_m=np.column_stack([xi[:-1],xi[1:]]),
                wall=mesh,cell_parent=parent,cell_branch=branch,points=points,sample_cell=sample_cell,
                sample_parent=parent[sample_cell],normals=np.column_stack([-tangent[sample_cell,1],tangent[sample_cell,0]]),
                q=q,cell_power=power,areas=2*np.pi*(mesh[:-1,0]+mesh[1:,0])/2*ds,
                ledger=ledger,lengths_m=L,sigma_m=sigma,source_xi_m=source_xi,valid=not reasons,reasons=reasons,
                diagnostics=dict(path_flux_drift_xi_m=drift,donor_branch=donor_branch,
                                 affected_blocked_tubes=int(np.count_nonzero(blocked>0)),
                                 affected_branch_tail_tubes=int(np.count_nonzero(tail>1e-14)),
                                 unhit_tubes=int(np.count_nonzero(active&~hit)),visibility_fans=visibility_tests,visibility_method="whole_cell_triangular_fan",
                                 gaussian_support_sigma=8.,config=asdict(config)))
