"""Independent first hits/thermal quadrature and fixed-candidate convergence."""
import argparse,time
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.integrate import quad
from shapely.geometry import LineString,Polygon
from d3d_optimization.boundary_evaluator import prepare_boundary,BoundaryEvaluator,export_boundary
from d3d_optimization.boundary_transport import AreaBoundary,wrapped_density
from d3d_optimization.reference_transport import ReferencePoint
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.broadening_study import load_json,save,digest
from verification.validate_reference_transport import independent_thermal
from optimize_boundary_transport import verify

LEVELS=[('cubic',320,.0005,161,.001),('cubic',640,.00025,321,.0005),
        ('cubic',1280,.000125,641,.00025),('quintic',640,.00025,321,.0005)]


def audit(e,r):
    row=e.describe(r);row['audit_passed']=False
    if 'deposition' not in r:return row
    d=r['deposition'];full=first_hits(e.paths,r['wall']);good=full.segment>=0
    mismatches=int(np.count_nonzero(full.segment!=r['hits'].segment))
    point_error=float(np.max(np.linalg.norm(full.points[good]-r['hits'].points[good],axis=1))) if np.any(good) else 0.
    independent=independent_thermal(e.field,e.paths,full,e.metadata['source_RZ_m'],e.thermal_paths.qu,e.reference,order=256)
    errors={k:float(np.max(abs(independent[k][good]-r['thermal'][k][good])/independent[k][good])) for k in independent}
    wall=LineString(r['wall'])
    geometry=bool(wall.is_simple and Polygon(r['wall']).is_valid and not wall.intersects(Polygon(e.lcfs))
                  and not wall.intersects(e.vessel_solids))
    fixed_error=0.
    for i in e.fixed_indices:
        ids=np.flatnonzero(r['original_segment']==i)
        fixed_error=max(fixed_error,float(np.max(abs(r['wall'][ids[0]:ids[0]+2]-e.wall[i:i+2]))))
    rsep,z=e.metadata['source_RZ_m'];sign=-1 if e.metadata['source']=='IMP' else 1
    def density(x):
        pt=[[rsep+sign*x,z]]
        return float(np.exp(-x/e.source.width)*2*np.pi*abs(e.field.gradient(pt)[0,0])/np.linalg.norm(e.field.magnetic(pt)[0]))
    integral=quad(density,0,e.source.edges[-1],epsabs=1e-13,epsrel=1e-10)[0]
    normalization_error=abs(integral-e.source.normalization)/integral
    # Independent dense donor-adaptive samples check the continuous peak search.
    mesh=AreaBoundary(r['wall'],e.wall_step)
    centers=d['donor_area_m2'][good];widths=d['sigma_area_m2'][good];weights=e.source.weights[good]
    area=np.r_[mesh.edges,(centers[:,None]+np.minimum(widths[:,None],mesh.total)*np.linspace(-3,3,49)).ravel()]%mesh.total
    q=np.zeros(len(area))
    for c,s,w in zip(centers,widths,weights):q+=w*wrapped_density(area,c,s,mesh.total)
    sample_peak=float(q.max());peak_underestimate=max(0.,sample_peak/r['peak_q']-1)
    arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(d['wall'],axis=0),axis=1))]
    cdf=np.interp(np.linspace(0,1,201),arc/arc[-1],np.r_[0.,np.cumsum(d['cell_power'])])
    row.update(full_first_hit_mismatches=mismatches,full_first_hit_point_error_m=point_error,
        full_wall_geometry_passed=geometry,fixed_wall_error_m=fixed_error,source_normalization_relative_error=normalization_error,
        independent_thermal_max_relative_errors=errors,independent_dense_peak_underestimate=peak_underestimate,
        independent_dense_peak_MW_m2=sample_peak*e.reference.power_W/1e6,cdf_wall_fraction=cdf.tolist())
    row['audit_passed']=bool(r['feasible'] and geometry and fixed_error<1e-12 and mismatches==0 and point_error<1e-10
        and normalization_error<1e-8 and max(errors.values())<.005 and peak_underestimate<1e-5
        and row['area_quadrature_error']<1e-4 and d['ledger']['balance_error']<2e-12
        and row['thermal']['conduction_relative_residual']<1e-10 and row['thermal']['sheath_relative_residual']<1e-10)
    return row


def verdict(levels,wall_only=None):
    result=dict(improvement_confirmed=False,resolution_screen_passed=False,global_optimum_proven=False,
                scope='same boundary-footprint model and fixed operating inputs; numerical comparison only')
    if len(levels)!=4 or not all(all(l.get('rows',{}).get(n,{}).get('audit_passed') and l['rows'][n].get('feasible')
        for n in ['baseline','selected']) for l in levels):return result
    if (wall_only is None or wall_only.get('max_peak_relative_change',float('inf'))>.005 or
        not all(wall_only.get('rows',{}).get(n,{}).get('audit_passed') and wall_only['rows'][n].get('feasible')
                for n in ['baseline','selected'])):return result
    ratios=np.array([l['rows']['selected']['peak_ratio'] for l in levels])
    peaks=np.array([[l['rows'][n]['peak_MW_m2'] for n in ['baseline','selected']] for l in levels])
    if not np.all(np.isfinite(peaks)) or not np.all(np.isfinite(ratios)) or np.any(peaks<=0):return result
    ratio_change=float(max(abs(ratios[2]-ratios[1]),abs(ratios[3]-ratios[1])))
    peak_change=float(max(np.max(abs(peaks[j]-peaks[1])/peaks[1]) for j in [2,3]))
    region_change=cdf_change=location_change=0.
    for n in ['baseline','selected']:
        base=levels[1]['rows'][n]
        for j in [2,3]:
            other=levels[j]['rows'][n]
            region_change=max(region_change,sum(abs(base['region_power'][k]-other['region_power'][k]) for k in base['region_power']))
            cdf_change=max(cdf_change,float(np.max(abs(np.array(base['cdf_wall_fraction'])-other['cdf_wall_fraction']))))
            location_change=max(location_change,float(np.linalg.norm(np.array(base['peak_RZ_m'])-other['peak_RZ_m'])))
    passed=ratio_change<=.002 and peak_change<=.005 and region_change<=.005 and cdf_change<=.002 and location_change<=.002
    gain=1-float(ratios[2])
    result.update(resolution_screen_passed=bool(passed),max_peak_ratio_change=ratio_change,
        max_absolute_peak_relative_change=peak_change,max_region_power_l1_change=region_change,
        max_wall_cdf_change=cdf_change,max_peak_location_change_m=location_change,finest_reduction_percent=100*gain,
        improvement_confirmed=bool(passed and np.all(ratios<=1+1e-10) and gain>2*ratio_change+1e-4))
    return result


def validate_source(folder):
    folder=Path(folder);verify(folder);manifest=load_json(folder/'manifest.json')
    if manifest.get('schema')!='boundary-transport-study.v1' or 'selected' not in manifest['candidates']:
        raise ValueError('selected boundary-study candidate required')
    if (folder/'validation.json').exists():raise FileExistsError('validation evidence already exists')
    record=dict(complete=False,levels=[],validator_sha256=digest(Path(__file__)),thresholds=dict(
        peak_relative=.005,ratio_absolute=.002,region_l1=.005,wall_cdf=.002,peak_location_m=.002))
    reference=ReferencePoint(**manifest['reference'])
    for index,(method,count,step,nodes,wall_step) in enumerate(LEVELS):
        start=time.perf_counter()
        e=prepare_boundary(ROOT,manifest['source'],reference,count=count,step=step,nodes=nodes,wall_step=wall_step,
            method=method,quadrature=manifest['settings']['quadrature'],bound=manifest['settings']['bound_mm']/1000)
        if any(e.metadata[k]!=manifest['metadata'][k] for k in ['shot_sha256','geometry_sha256','surfaces']):raise ValueError('validation input mismatch')
        rows={}
        for name in ['baseline','selected']:
            r=e.baseline if name=='baseline' else e.evaluate(manifest['candidates'][name])
            rows[name]=audit(e,r)
            if index==2:export_boundary(folder/(name+'_refined.npz'),e,r)
        level=dict(method=method,tubes=count,trace_step_m=step,nodes=nodes,wall_step_m=wall_step,rows=rows,
                   elapsed_seconds=time.perf_counter()-start)
        record['levels'].append(level);save(folder/'validation.json',record)
        print(manifest['source'],method,count,'audit',[r['audit_passed'] for r in rows.values()],
              'peaks',[r.get('peak_MW_m2') for r in rows.values()],flush=True)
        if index==0:
            fine=BoundaryEvaluator(e.direct,e.metadata,reference,wall_step/2,thermal_paths=e.thermal_paths)
            wall_rows={n:audit(fine,fine.baseline if n=='baseline' else fine.evaluate(manifest['candidates'][n])) for n in rows}
            record['wall_only_refinement']=dict(wall_step_m=wall_step/2,rows=wall_rows,
                max_peak_relative_change=max(abs(wall_rows[n]['peak_q']/rows[n]['peak_q']-1) for n in rows))
    record['verdict']=verdict(record['levels'],record['wall_only_refinement'])
    # Require the production search screen and independent receiving-grid check too.
    if (manifest['settings']['smoke'] or not manifest['baseline']['feasible'] or not manifest['selected']['feasible']
            or record['wall_only_refinement']['max_peak_relative_change']>.005):
        record['verdict']['improvement_confirmed']=False
    record['complete']=True;save(folder/'validation.json',record)
    manifest['validation_status']='complete_see_validation_verdict';save(folder/'manifest.json',manifest)
    return record


def validate_study(output,source='OMP'):
    return {s:validate_source(Path(output)/s)['verdict'] for s in ['IMP','OMP'] if source=='both' or s==source}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--source',default='OMP',choices=['IMP','OMP','both'])
    a=p.parse_args();print(validate_study(a.output,a.source))
