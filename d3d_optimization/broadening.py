"""Power-accounted connection-length broadening in upstream-equivalent flux.

This phenomenological endpoint surrogate is not a perpendicular transport solver.
Gaussian widths are standard deviations. No lost/blocked power is renormalized.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.special import ndtr
from scipy.optimize import brentq
import shapely
from shapely.geometry import LineString, Polygon


@dataclass(frozen=True)
class BroadeningConfig:
    sigma_ref_m: float
    length_ref_m: float = 10.
    exponent: float = 1.
    wall_step_m: float = .002
    anchor_m: float = .005
    missing_tolerance: float = .001

    def __post_init__(self):
        if (not all(np.isfinite(v) for v in asdict(self).values()) or
                self.sigma_ref_m < 0 or self.length_ref_m <= 0 or
                self.exponent not in (.5, 1.) or self.wall_step_m <= 0 or
                self.anchor_m <= 0 or not 0 <= self.missing_tolerance < 1):
            raise ValueError('invalid finite broadening parameters')


class PathLengths:
    """Cached trapezoidal integral of |B|/|Bp| on actual poloidal paths.

    A fractional segment uses the integral of its linearly interpolated metric,
    rather than a linear interpolation of cumulative length. Invalid metrics
    fail on use; metrics beyond a first hit cannot invalidate that first hit.
    """
    def __init__(self, field, paths):
        self.paths = [np.asarray(p.coords if hasattr(p, 'coords') else p, float) for p in paths]
        self.s, self.metric, self.cumulative, self.psi = [], [], [], []
        for p in self.paths:
            if p.ndim != 2 or p.shape[1] != 2 or len(p)<2 or not np.all(np.isfinite(p)):
                raise ValueError('finite paths with at least two vertices required')
            ds=np.linalg.norm(np.diff(p,axis=0),axis=1)
            if np.any(ds<=0): raise ValueError('duplicate path vertices forbidden')
            b=np.asarray(field.magnetic(p),float)
            bp=np.linalg.norm(b[:,:2],axis=1)
            with np.errstate(invalid='ignore',divide='ignore'):
                metric=np.linalg.norm(b,axis=1)/bp
            self.s.append(np.r_[0.,np.cumsum(ds)])
            self.metric.append(metric)
            self.cumulative.append(np.r_[0.,np.cumsum(ds*(metric[:-1]+metric[1:])/2)])
            self.psi.append(np.asarray(field.values(p),float))

    def at(self, distances):
        d=np.asarray(distances,float)
        if d.shape != (len(self.paths),): raise ValueError('one distance per path required')
        result=np.full(len(d),np.nan)
        for i,x in enumerate(d):
            if np.isposinf(x): continue
            s=self.s[i]
            if not np.isfinite(x) or x<0 or x>s[-1]+1e-10:
                raise ValueError('finite hit must lie on cached path')
            if x==0:
                result[i]=0.
                continue
            x=min(x,s[-1]); j=min(np.searchsorted(s,x,side='left')-1,len(s)-2)
            m=self.metric[i]
            if not np.all(np.isfinite(m[:j+2])) or np.any(m[:j+2]<1):
                raise ValueError('invalid Bp or nonfinite connection length')
            h=x-s[j]; f=h/(s[j+1]-s[j])
            result[i]=self.cumulative[i][j]+h*(m[j]+.5*f*(m[j+1]-m[j]))
        return result


def _cdf_interval(lo, hi):
    """Stable normal probability, including intervals in the positive tail."""
    return np.where(lo>=0, ndtr(-lo)-ndtr(-hi), ndtr(hi)-ndtr(lo))


def _wall_mesh(field, wall, scale, offset, step):
    wall=np.asarray(wall,float)
    if (wall.ndim!=2 or wall.shape[1]!=2 or len(wall)<4 or
            not np.all(np.isfinite(wall)) or np.any(wall[:,0]<=0) or
            np.linalg.norm(wall[0]-wall[-1])>1e-10):
        raise ValueError('finite closed positive-R wall required')
    # On each tensor-spline patch the directional derivative along a straight
    # segment is a polynomial of degree kr+kz-1 (5 cubic, 9 quintic). Split at
    # every grid/knots crossing before reconstructing that polynomial, so two
    # close extrema cannot hide between finite sign probes.
    if getattr(field,'method',None)=='linear':
        knot_r,knot_z=field.r,field.z
        degree=1
    elif hasattr(field,'spline'):
        knot_z,knot_r=field.spline.get_knots()
        degree=sum(field.spline.degrees)-1
    else:
        knot_r,knot_z=[],[]
        degree=9  # Generic test fields must pass the polynomial residual check.
    points=[wall[0]]; parents=[]
    for parent,(a,b) in enumerate(zip(wall[:-1],wall[1:])):
        length=np.linalg.norm(b-a)
        if length<=0: raise ValueError('duplicate wall vertices forbidden')
        cuts=[np.linspace(0,1,int(np.ceil(length/step))+1)]
        for axis,knots in enumerate((knot_r,knot_z)):
            if b[axis]!=a[axis]:
                f=(np.unique(knots)-a[axis])/(b[axis]-a[axis])
                cuts.append(f[(f>0)&(f<1)])
        fractions=np.unique(np.concatenate(cuts))
        fractions=fractions[np.r_[True,np.diff(fractions)>1e-13]][1:]
        if fractions[-1]<1: fractions=np.r_[fractions,1.]
        points.extend(a+fractions[:,None]*(b-a))
        parents.extend([parent]*len(fractions))
    initial=np.asarray(points)
    delta=np.diff(initial,axis=0)
    tangent0=delta/np.linalg.norm(delta,axis=1)[:,None]
    # Interior Chebyshev roots avoid evaluating a discontinuous bilinear
    # derivative on the wrong side of a patch endpoint.
    fit_nodes=np.cos(np.pi*(np.arange(degree+1)+.5)/(degree+1))
    check_nodes=np.cos(np.pi*(np.arange(degree+2)+.5)/(degree+2))
    nodes=np.r_[fit_nodes,check_nodes]
    fractions=(nodes+1)/2
    probe=initial[:-1,None,:]+fractions[None,:,None]*delta[:,None,:]
    gradients=scale*np.asarray(field.gradient(probe.reshape(-1,2))).reshape(-1,len(nodes),2)
    derivatives=(gradients*tangent0[:,None,:]).sum(axis=2)
    if not np.all(np.isfinite(derivatives)): raise ValueError('nonfinite receiving gradient')
    vandermonde=np.polynomial.chebyshev.chebvander(fit_nodes,degree)
    coefficients=derivatives[:,:degree+1]@np.linalg.inv(vandermonde).T
    predicted=coefficients@np.polynomial.chebyshev.chebvander(check_nodes,degree).T
    # Include component magnitudes in the floating-point allowance when the
    # directional projection cancels near a tangent wall.
    component_scale=np.max(np.abs(gradients),axis=(1,2))
    local_scale=np.max(np.abs(derivatives),axis=1)
    tolerance=1e-10*local_scale+256*np.finfo(float).eps*component_scale+1e-15
    if np.any(np.max(abs(predicted-derivatives[:,degree+1:]),axis=1)>tolerance):
        raise ValueError('receiving derivative is not a resolved patch polynomial')
    # The Chebyshev bound certifies a fixed sign for most cells. Only the
    # remainder need a companion-matrix root solve; all roots are checked.
    bounded=abs(coefficients[:,0])>np.sum(abs(coefficients[:,1:]),axis=1)+tolerance
    roots_by_cell={}
    for j in np.flatnonzero(~bounded):
        c=coefficients[j].copy()
        cutoff=64*np.finfo(float).eps*max(float(np.max(abs(c))),1e-300)
        significant=np.flatnonzero(abs(c)>cutoff)
        if len(significant)<1 or significant[-1]==0: continue
        roots=np.polynomial.chebyshev.chebroots(c[:significant[-1]+1])
        real=roots.real[(abs(roots.imag)<1e-7)&(roots.real>-1+1e-10)&(roots.real<1-1e-10)]
        f=np.sort((real+1)/2)
        f=f[np.r_[True,np.diff(f)>1e-10]] if len(f) else f
        if len(f):
            p=initial[j]+f[:,None]*delta[j]
            residual=abs(scale*np.asarray(field.gradient(p))@tangent0[j])
            if np.any(residual>max(1e-8*local_scale[j],10*tolerance[j])):
                raise ValueError('unresolved receiving polynomial root')
            roots_by_cell[j]=f
    if roots_by_cell:
        new_points=[initial[0]]; new_parents=[]
        for j,(a,b) in enumerate(zip(initial[:-1],initial[1:])):
            for end in np.r_[roots_by_cell.get(j,[]),1.]:
                new_points.append(a+end*(b-a)); new_parents.append(parents[j])
        points,parents=new_points,new_parents
    points=np.asarray(points); parent=np.asarray(parents,int)
    d=np.diff(points,axis=0); ds=np.linalg.norm(d,axis=1); tangent=d/ds[:,None]
    xi=scale*(np.asarray(field.values(points))-offset)
    probe=points[:-1,None,:]+np.array([.1,.3,.5,.7,.9])[None,:,None]*d[:,None,:]
    derivatives=(scale*np.asarray(field.gradient(probe.reshape(-1,2))).reshape(-1,5,2)*tangent[:,None,:]).sum(axis=2)
    if not np.all(np.isfinite(xi)) or not np.all(np.isfinite(derivatives)):
        raise ValueError('nonfinite flux mesh')
    tol=1e-12*max(1.,float(np.max(np.abs(derivatives))))
    if np.any((np.min(derivatives,axis=1)<-tol)&(np.max(derivatives,axis=1)>tol)):
        raise ValueError('unresolved wall flux extremum; refine receiving grid')
    sign=np.sign(np.diff(xi)).astype(int)
    if np.any((np.diff(xi)*np.mean(derivatives,axis=1)) < -tol):
        raise ValueError('inconsistent wall flux monotonicity')
    branch=np.cumsum(np.r_[0,(sign[1:]!=sign[:-1]) | (sign[1:]==0)])
    return points,parent,ds,tangent,xi,branch


def build_deposition(field, paths, lengths, hits, wall, weights_matrix,
                     source_rz, source_sign, lcfs, vessel_solids, config):
    """Deposit all source widths in a single geometry/kernel pass.

    `source_rz` is the separatrix point. Actual donor flux is sampled at each
    path start. Weights have shape (donor, source width). Parent segment IDs
    refer to supplied candidate wall segments, not original design segments.
    """
    if config.sigma_ref_m<=0: raise ValueError('zero spread must dispatch to direct evaluator')
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
    L=lengths.at(distance); sigma=config.sigma_ref_m*(L/config.length_ref_m)**config.exponent
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
