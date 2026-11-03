"""CPD-5 free-boundary equilibrium and immutable sampled field.

All controlled currents are total ampere-turns; electrical turns are unknown.
The pressure fit preserves the explicitly declared reference beta_N convention.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import yaml
from scipy.interpolate import RectBivariateSpline
import freegs
from freegs.machine import Machine, Solenoid
from freegs.shaped_coil import ShapedCoil

MU0=4e-7*np.pi

def reference_inputs():
    return dict(Ip_A=1388124.0684837103,Bt_T=2e-7*144*123784.60159301758/1.715993046760559,
                betaN=0.6166043877601624,R0_m=4.007,a_reference_m=1.272,
                alpha_m=1.,alpha_n=2.,D3D_time_ms=1400.,D3D_bcoil_kA=123.78460159301758,
                D3D_Raxis_m=1.715993046760559,D3D_tf_turns=144,
                betaN_definition='100 * (2 mu0 volume-average pressure / Bt_reference_T^2) * measured_a_LCFS * Bt_reference_T / Ip_MA; approximate D3D operating match, EFIT convention not fully known')

def normalized_beta(pressure_average, minor_radius, inputs):
    return 100.*(2*MU0*pressure_average/inputs['Bt_T']**2)*minor_radius*inputs['Bt_T']/(inputs['Ip_A']/1e6)

class DistributedSolenoid(Solenoid):
    """Uniform axial filament quadrature normalized to one total ampere-turn.

    Ns is numerical quadrature size, never an assumed electrical turn count.
    The physical winding radius is the mean of the YAML OH pack edges.
    """
    def __init__(self, radius,zmin,zmax,quadrature_points=101):
        super().__init__(radius,zmin,zmax,quadrature_points,current=0.,control=True)
        self.turns=1
    def controlPsi(self,R,Z): return super().controlPsi(R,Z)/self.Ns
    def controlBr(self,R,Z): return super().controlBr(R,Z)/self.Ns
    def controlBz(self,R,Z): return super().controlBz(R,Z)/self.Ns

class FrozenField:
    """Read NPZ psi[Z,R] once; no equilibrium or coil updates during tracing."""
    def __init__(self,path):
        self.path=Path(path)
        self.sha256=hashlib.sha256(self.path.read_bytes()).hexdigest()
        with np.load(path,allow_pickle=False) as d:
            self.R=d['R'].copy(); self.Z=d['Z'].copy(); self.psi=d['psi'].copy()
            self.fvac=float(d['fvac']); self.psi_axis=float(d['psi_axis']); self.psi_lcfs=float(d['psi_lcfs'])
        if self.psi.shape != (len(self.Z),len(self.R)) or not np.isfinite(self.psi).all():
            raise ValueError('Expected finite psi[Z,R]')
        if not (np.all(np.diff(self.R)>0) and np.all(np.diff(self.Z)>0)):
            raise ValueError('Grid axes must increase strictly')
        self._spline=RectBivariateSpline(self.R,self.Z,self.psi.T)
        self.r,self.z=self.R,self.Z
        for x in (self.R,self.Z,self.psi): x.setflags(write=False)
    def checked(self,points):
        p=np.asarray(points,dtype=float).reshape(-1,2)
        self._check(p[:,0],p[:,1]);return p
    def values(self,points):
        p=self.checked(points);return self.flux(p[:,0],p[:,1])
    def gradient(self,points):
        p=self.checked(points)
        return np.column_stack([self._spline.ev(p[:,0],p[:,1],dx=1),self._spline.ev(p[:,0],p[:,1],dy=1)])
    def magnetic(self,points):
        p=self.checked(points);return np.column_stack(self.field(p[:,0],p[:,1]))
    def incidence(self,points,normals):
        b=self.magnetic(points);n=np.asarray(normals,float)
        if n.shape!=b[:,:2].shape or not np.isfinite(n).all() or np.any(np.linalg.norm(n,axis=1)==0):
            raise ValueError('Finite nonzero R/Z normals required')
        n=n/np.linalg.norm(n,axis=1)[:,None]
        return np.degrees(np.arcsin(np.clip(abs(np.sum(b[:,:2]*n,axis=1))/np.linalg.norm(b,axis=1),0,1)))
    def _check(self,r,z):
        r,z=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(z,dtype=float))
        if not (np.isfinite(r).all() and np.isfinite(z).all()) or np.any((r<self.R[0])|(r>self.R[-1])|(z<self.Z[0])|(z>self.Z[-1])):
            raise ValueError('Point outside frozen field domain')
        return r,z
    def flux(self,r,z):
        r,z=self._check(r,z);return self._spline.ev(r,z)
    def field(self,r,z):
        r,z=self._check(r,z)
        return -self._spline.ev(r,z,dy=1)/r,self._spline.ev(r,z,dx=1)/r,self.fvac/r

def build_machine(nodes):
    data=yaml.safe_load(Path(nodes).read_text(encoding='utf-8'))
    if data['metadata']['units']!='mm': raise ValueError('CPD source must declare mm')
    coils=[];details={}
    for name,spec in data['components'].items():
        if not (name.startswith('PF_') or name=='OH_Solenoid'): continue
        p=np.array([[n['x'],n['z']] for n in spec['contours']['outer']['nodes']],float)*.001
        if name.startswith('PF_'):
            coil=ShapedCoil(p.tolist(),turns=1,control=True,npoints=6)
            model='FreeGS ShapedCoil, six quadrature points per triangle'
        else:
            coil=DistributedSolenoid(float(p[:,0].mean()),float(p[:,1].min()),float(p[:,1].max()))
            model='101 axial filaments normalized to total ampere-turns at mean OH radius'
        coils.append((name,coil))
        details[name]=dict(outline_m=p.tolist(),model=model,effective_turns=1,electrical_turns=None)
    return Machine(coils,R0=4.007),details

def target_shape():
    return dict(xpoint_m=[3.55,-1.80],isoflux_points_m=[[2.90,0.],[5.10,0.],[3.55,1.60],[4.70,1.10]],
                domain_m=[1.65,7.60,-4.10,4.10],
                rationale='Geometry-derived lower single null; 0.165/0.179m midplane blanket gaps; upper top (3.55,1.60) and shoulder (4.70,1.10) chosen after taller diagnostic targets produced upper-limiting topology; no input Plasma contour used.')

def seed_equilibrium(machine,inputs,target,n):
    rmin,rmax,zmin,zmax=target['domain_m']
    r,z=np.meshgrid(np.linspace(rmin,rmax,n),np.linspace(zmin,zmax,n),indexing='ij')
    for _,coil in machine.coils: coil.current=0.
    psi=MU0*inputs['Ip_A']*inputs['R0_m']/4*np.exp(-((r-inputs['R0_m'])/(1.3*inputs['a_reference_m']))**2-((z-.1)/(2*inputs['a_reference_m']))**2)
    psi[0,:]=psi[-1,:]=psi[:,0]=psi[:,-1]=0.
    return freegs.Equilibrium(tokamak=machine,Rmin=rmin,Rmax=rmax,Zmin=zmin,Zmax=zmax,nx=n,ny=n,psi=psi,boundary=freegs.boundary.freeBoundaryHagenow)

def solve_pressure(nodes,n=65,paxis=20000.,inputs=None,target=None,eq=None,history=None):
    inputs=reference_inputs() if inputs is None else inputs
    target=target_shape() if target is None else target
    history=[] if history is None else history
    if eq is None:
        machine,_=build_machine(nodes);eq=seed_equilibrium(machine,inputs,target,n)
        stages=[(None,1e-7,.65)]+[(target['xpoint_m'],g,b) for g,b in zip((1e-7,3e-8,1e-8,3e-9,1e-9),(.65,.55,.45,.35,.25))]
    else: stages=[(target['xpoint_m'],1e-9,.35)]
    profiles=freegs.jtor.ConstrainPaxisIp(eq,paxis,inputs['Ip_A'],inputs['R0_m']*inputs['Bt_T'],alpha_m=inputs['alpha_m'],alpha_n=inputs['alpha_n'],Raxis=1.)
    iso=[tuple(target['xpoint_m']+p) for p in target['isoflux_points_m']]
    for xp,gamma,blend in stages:
        control=freegs.control.constrain(xpoints=[] if xp is None else [xp],isoflux=iso,gamma=gamma)
        rtol=2e-3 if gamma>1e-9 else 3e-5
        convergence=freegs.solve(eq,profiles,control,rtol=rtol,blend=blend,maxits=400,convergenceInfo=True)
        record=dict(n=n,paxis_Pa=paxis,gamma=gamma,blend=blend,Ip_A=float(eq.plasmaCurrent()),iterations=len(convergence[1]),psi_relative_change=float(convergence[1][-1]))
        history.append(record)
        print(json.dumps(record),flush=True)
    return eq,history

def measure(eq,inputs,target,geometry=None):
    op,xp=freegs.critical.find_critical(eq.R,eq.Z,eq.psi())
    if not op or not xp: raise RuntimeError('No magnetic axis or X-point')
    lcfs=np.asarray(eq.separatrix(npoints=720))
    from shapely.geometry import Polygon
    poly=Polygon(lcfs)
    beta=normalized_beta(float(eq.pressure_ave()),float(eq.minorRadius()),inputs)
    result=dict(Ip_A=float(eq.plasmaCurrent()),Ip_relative_error=float(abs(eq.plasmaCurrent()/inputs['Ip_A']-1)),betaN=beta,betaN_relative_error=float(abs(beta/inputs['betaN']-1)),betaN_FreeGS=float(eq.betaN()),paxis_Pa=float(eq.pressure(0.)),axis_m=list(map(float,op[0][:2])),xpoints=[list(map(float,p)) for p in xp],psi_axis=float(eq.psi_axis),psi_lcfs=float(eq.psi_bndry),lower_single_null=bool(xp[0][1]<0 and all(abs((p[2]-xp[0][2])/(eq.psi_bndry-eq.psi_axis))>.01 for p in xp[1:])),lcfs_valid_polygon=bool(poly.is_valid),lcfs_bounds=list(poly.bounds))
    result['minor_radius_m']=float(eq.minorRadius())
    result['pressure_average_Pa']=float(eq.pressure_ave())
    result['xpoint_target_distance_m']=float(np.linalg.norm(np.array(xp[0][:2])-target['xpoint_m']))
    if geometry is not None:
        result['material_overlap_area_m2']={k:float(poly.intersection(s).area) for k,s in geometry.solids.items()}
        result['material_clearance_m']={k:float(poly.distance(s)) for k,s in geometry.solids.items()}
    return result,lcfs

def save_solution(eq,nodes,output,inputs,target,history,geometry=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    metrics,lcfs=measure(eq,inputs,target,geometry)
    _,coils=build_machine(nodes)
    for name,coil in eq.tokamak.coils: coils[name]['total_ampere_turns']=float(coil.current)
    np.savez_compressed(output/'equilibrium.npz',R=eq.R[:,0],Z=eq.Z[0,:],psi=eq.psi().T,fvac=inputs['R0_m']*inputs['Bt_T'],psi_axis=eq.psi_axis,psi_lcfs=eq.psi_bndry,lcfs=lcfs,axis=metrics['axis_m'],xpoint=metrics['xpoints'][0][:2])
    frozen=FrozenField(output/'equilibrium.npz')
    sr,sz=np.meshgrid(np.linspace(2.7,5.4,25),np.linspace(-2.6,2.6,31))
    points=np.column_stack([sr.ravel(),sz.ravel()])
    direct=np.column_stack([eq.Br(points[:,0],points[:,1]),eq.Bz(points[:,0],points[:,1])])
    sampled=frozen.magnetic(points)[:,:2]
    error=np.linalg.norm(sampled-direct,axis=1)
    keep=np.linalg.norm(direct,axis=1)>.01
    metrics['frozen_field_parity']=dict(samples=len(points),max_Bp_absolute_error_T=float(error.max()),p95_Bp_relative_error_away_from_null=float(np.quantile(error[keep]/np.linalg.norm(direct[keep],axis=1),.95)))
    report=dict(status='diagnostic_not_validated',inputs=inputs,target=target,metrics=metrics,coils=coils,history=history,grid_order='psi[Z,R]',field_convention='Br=-dpsi/dZ/R; Bz=dpsi/dR/R; SOL Bphi=fvac/R',source_sha256=hashlib.sha256(Path(nodes).read_bytes()).hexdigest(),freegs_source=freegs.__file__,frozen_field_sha256=hashlib.sha256((output/'equilibrium.npz').read_bytes()).hexdigest())
    (output/'equilibrium.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,9))
    if geometry is not None:
        for name,s in geometry.solids.items():
            for part in (list(s.geoms) if hasattr(s,'geoms') else [s]):
                from shapely.geometry.polygon import orient
                part=orient(part,sign=1.)
                from matplotlib.path import Path as MplPath
                from matplotlib.patches import PathPatch
                rings=[MplPath(np.asarray(part.exterior.coords))]+[MplPath(np.asarray(r.coords)) for r in part.interiors]
                color=dict(Blanket='tab:blue',Divertor='tab:orange',Shield='tab:green',VV='tab:red')[name]
                ax.add_patch(PathPatch(MplPath.make_compound_path(*rings),facecolor=color,edgecolor=color,alpha=.25,label=name))
    pn=(eq.psi()-eq.psi_axis)/(eq.psi_bndry-eq.psi_axis)
    ax.contour(eq.R,eq.Z,pn,levels=[.2,.4,.6,.8,1,1.01,1.03,1.06],linewidths=.6)
    ax.plot(lcfs[:,0],lcfs[:,1],'k-',lw=1.4,label='LCFS')
    ax.scatter(*metrics['axis_m'],marker='+',c='red');ax.scatter(*metrics['xpoints'][0][:2],marker='x',c='red')
    ax.set(xlabel='R [m]',ylabel='Z [m]',title=f"CPD-5 diagnostic: Ip={metrics['Ip_A']/1e6:.4f} MA, betaN={metrics['betaN']:.5f}")
    ax.set_aspect('equal');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(output/'equilibrium.png',dpi=160);plt.close(fig)
    return report








def equilibrium_failures(metrics):
    """Numerical prerequisite checks; no engineering current limit is assumed."""
    failures=[]
    for name,key,limit in [('Ip','Ip_relative_error',1e-5),('betaN','betaN_relative_error',.002),('xpoint_control','xpoint_target_distance_m',.005)]:
        value=metrics.get(key,float('inf'))
        if not np.isfinite(value) or value>limit: failures.append(name)
    for key in ('lower_single_null','lcfs_valid_polygon'):
        if not metrics.get(key,False): failures.append(key)
    if not metrics.get('material_overlap_area_m2') or not metrics.get('material_clearance_m'):
        failures.append('material_checks_missing')
    else:
        if any(not np.isfinite(v) or v>1e-10 for v in metrics['material_overlap_area_m2'].values()): failures.append('material_overlap')
        if any(not np.isfinite(v) or v<.005 for v in metrics['material_clearance_m'].values()): failures.append('material_clearance')
    parity=metrics.get('frozen_field_parity',{}).get('p95_Bp_relative_error_away_from_null',float('inf'))
    if not np.isfinite(parity) or parity>.005: failures.append('frozen_field_parity')
    return failures

def compare_refinement(coarse_dir,fine_dir):
    """Compare independently converged target solves on nested uniform grids."""
    from shapely.geometry import Polygon
    coarse_dir,fine_dir=Path(coarse_dir),Path(fine_dir)
    c=json.loads((coarse_dir/'equilibrium.json').read_text());f=json.loads((fine_dir/'equilibrium.json').read_text())
    cf,ff=FrozenField(coarse_dir/'equilibrium.npz'),FrozenField(fine_dir/'equilibrium.npz')
    failures=equilibrium_failures(c['metrics'])+equilibrium_failures(f['metrics'])
    if len(ff.r)<=len(cf.r) or len(ff.z)<=len(cf.z): failures.append('grid_not_refined')
    if c['inputs']!=f['inputs'] or c['target']!=f['target'] or c['source_sha256']!=f['source_sha256']: failures.append('inputs_changed')
    r,z=np.meshgrid(np.linspace(2.7,5.4,151),np.linspace(-2.2,2.,191))
    p=np.column_stack([r.ravel(),z.ravel()])
    cn=(cf.values(p)-cf.psi_axis)/(cf.psi_lcfs-cf.psi_axis)
    fn=(ff.values(p)-ff.psi_axis)/(ff.psi_lcfs-ff.psi_axis)
    cb,fb=cf.magnetic(p)[:,:2],ff.magnetic(p)[:,:2]
    magnitude=np.linalg.norm(fb,axis=1)
    select=(fn>=.95)&(fn<=1.1)&(magnitude>.01)
    if not select.any(): raise ValueError('No near-separatrix comparison samples')
    with np.load(coarse_dir/'equilibrium.npz') as a: cl=Polygon(a['lcfs'])
    with np.load(fine_dir/'equilibrium.npz') as a: fl=Polygon(a['lcfs'])
    delta=np.linalg.norm(cb[select]-fb[select],axis=1)
    ci=np.array([c['coils'][k]['total_ampere_turns'] for k in sorted(c['coils'])])
    fi=np.array([f['coils'][k]['total_ampere_turns'] for k in sorted(f['coils'])])
    metrics=dict(coarse_grid=[len(cf.r),len(cf.z)],fine_grid=[len(ff.r),len(ff.z)],near_separatrix_samples=int(select.sum()),Bp_p95_relative_difference=float(np.quantile(delta/magnitude[select],.95)),Bp_max_absolute_difference_T=float(delta.max()),psiN_max_difference=float(np.max(abs(cn[select]-fn[select]))),lcfs_hausdorff_m=float(cl.hausdorff_distance(fl)),lcfs_area_relative_difference=float(abs(cl.area/fl.area-1)),axis_displacement_m=float(np.linalg.norm(np.array(c['metrics']['axis_m'])-f['metrics']['axis_m'])),xpoint_displacement_m=float(np.linalg.norm(np.array(c['metrics']['xpoints'][0][:2])-f['metrics']['xpoints'][0][:2])),coil_max_change_relative_to_max_current=float(np.max(abs(ci-fi))/np.max(abs(fi))))
    limits=dict(Bp_p95_relative_difference=.01,Bp_max_absolute_difference_T=.002,psiN_max_difference=.002,lcfs_hausdorff_m=.002,lcfs_area_relative_difference=.01,axis_displacement_m=.002,xpoint_displacement_m=.001,coil_max_change_relative_to_max_current=.05)
    failures.extend(key for key,limit in limits.items() if not np.isfinite(metrics[key]) or metrics[key]>limit)
    result=dict(status='validated_for_frozen_transport' if not failures else 'failed_prerequisite',failures=failures,metrics=metrics,thresholds=limits,coarse_field_sha256=cf.sha256,frozen_field_sha256=ff.sha256,limitations=['No winding electrical turn counts or engineering current-density limits supplied; currents are total ampere-turns.','D3D betaN convention is approximated by the documented vacuum-reference definition.','Grid comparison supports this prescribed equilibrium and is not a proof of continuum accuracy.'])
    (fine_dir/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    f['status']=result['status'];f['refinement_validation']=result
    (fine_dir/'equilibrium.json').write_text(json.dumps(f,indent=2),encoding='utf-8')
    return result

def verify_reference_source(path):
    """Read the one-row D3D parquet and verify every transferred scalar."""
    import d3d_optimization  # Registers the repository's optional PyArrow package.
    import pyarrow.parquet as pq
    path=Path(path)
    columns=['source','magnetics_plasma_current','magnetics_plasma_current_times','magnetics_bcoil','magnetics_time','efit_times','efit_beta_n','efit_r_axis']
    rows=pq.read_table(path,columns=columns).to_pylist()
    if len(rows)!=1: raise ValueError('Expected one D3D shot per parquet')
    s=rows[0];t=1400.
    idx=np.flatnonzero(np.isclose(s['efit_times'],t,rtol=0,atol=1e-6))
    if len(idx)!=1: raise ValueError('Exactly one EFIT 1400ms frame required')
    i=int(idx[0])
    for key in ('magnetics_plasma_current_times','magnetics_time'):
        times=np.asarray(s[key])
        if not np.isfinite(times).all() or np.any(np.diff(times)<=0) or not times[0]<=t<=times[-1]:
            raise ValueError('Reference interpolation requires increasing times bracketing 1400ms')
    values=dict(Ip_A=1000.*float(np.interp(t,s['magnetics_plasma_current_times'],s['magnetics_plasma_current'])),D3D_bcoil_kA=float(np.interp(t,s['magnetics_time'],s['magnetics_bcoil'])),betaN=float(s['efit_beta_n'][i]),D3D_Raxis_m=float(s['efit_r_axis'][i]))
    expected=reference_inputs()
    for key,value in values.items():
        if not np.isfinite(value) or not np.isclose(value,expected[key],rtol=1e-12,atol=1e-12):
            raise ValueError(f'D3D source value changed: {key}')
    return dict(path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source=s['source'],time_ms=t,efit_index=i,columns=columns,values=values,values_verified=True,magnetics_sampling='Linear interpolation at 1400ms on each signal time base',EFIT_sampling='Exact frame index 60 at 1400ms',units=dict(magnetics_plasma_current='kA',magnetics_bcoil='kA',times='ms',efit_r_axis='m',efit_beta_n='dimensionless normalized beta'),Bt_formula='2e-7 * 144 * bcoil_kA * 1000 / D3D_Raxis_m; same Bt at CPD reference R0=4.007m')
