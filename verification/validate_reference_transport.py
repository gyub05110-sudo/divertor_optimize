"""Independent magnetic/thermal quadrature and fixed-reference refinement."""
import argparse
from dataclasses import replace
from pathlib import Path
import sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.integrate import cumulative_simpson
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.independent_sources import prepare_evaluator
from d3d_optimization.reference_transport import ReferencePoint,solve_two_point,E_CHARGE
from d3d_optimization.transport_evaluator import TransportEvaluator,summarize
from d3d_optimization.reference_study import verify,acceptance,export
from d3d_optimization.broadening_study import load_json,save
from verification.validate_broadening import audit as geometry_audit

LEVELS=[('cubic',320,.0005,161,.001),('cubic',640,.00025,321,.0005),
        ('cubic',1280,.000125,641,.00025),('quintic',640,.00025,321,.0005)]


def independent_thermal(field,paths,hits,source_rz,qu,reference,order=512):
    """Re-sample actual cut paths, Simpson-integrate B, invert via dense samples.

    This never uses the production ThermalPaths cache or its partial-segment
    inversion. The physical two-point root is shared and checked analytically.
    """
    out={k:np.full(len(paths),np.nan) for k in ['sigma_m','length_m','integral_Bratio_m']}
    gradient=abs(float(field.gradient([source_rz])[0,0]))
    nodes,weights=np.polynomial.legendre.leggauss(order);u=(nodes+1)/2
    for i,(path,distance,segment) in enumerate(zip(paths,hits.distance,hits.segment)):
        if segment<0:continue
        p=np.asarray(path.coords);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
        cut=np.searchsorted(arc,distance,side='left')
        p=np.vstack([p[:cut],hits.points[i]])
        if len(p)<2:raise ValueError('independent positive path required')
        dense=np.empty((2*len(p)-1,2));dense[::2]=p;dense[1::2]=(p[:-1]+p[1:])/2
        ds=np.linalg.norm(np.diff(dense,axis=0),axis=1);s=np.r_[0,np.cumsum(ds)]
        magnetic=field.magnetic(dense);b=np.linalg.norm(magnetic,axis=1)
        metric=b/np.linalg.norm(magnetic[:,:2],axis=1)
        cumulative=cumulative_simpson(metric*b/b[0],x=s,initial=0)
        total=float(cumulative[-1])
        state=solve_two_point(float(qu[i]),b[-1]/b[0],total,reference)
        yu,yt=state['Tu_eV']**1.75,state['Tt_eV']**1.75
        a=3.5*qu[i]/reference.kappa0
        y=yt+u*a*total/(yu+yt)
        position=np.interp(total*(1-u*(y+yt)/(yu+yt)),cumulative,s)
        pts=np.column_stack([np.interp(position,s,dense[:,axis]) for axis in [0,1]])
        mapping=np.linalg.norm(field.gradient(pts),axis=1)/gradient
        bp=np.linalg.norm(field.magnetic(pts),axis=1)
        c=np.sqrt(1.5*E_CHARGE*reference.density_m3*state['Tu_eV']*reference.chi_perp_m2_s/reference.kappa0)
        out['sigma_m'][i]=2*c*total/(yu+yt)*np.dot(weights/2,mapping*b[0]/bp)/np.sqrt(2)
        out['integral_Bratio_m'][i]=total
        out['length_m'][i]=cumulative_simpson(metric,x=s,initial=0)[-1]
    return out


def audit(e,result,baseline_ok=True):
    # Existing independent wall/normalization/first-hit checks accept one width.
    row=geometry_audit(e,result,baseline_ok)
    row.update({k:v for k,v in summarize(e,result).items() if k in ['reference','mode','thermal','peak_q_MW_m2','deposited_width_statistics']})
    if 'thermal' not in result:return row
    hits=first_hits(e.paths,result['wall'])
    independent=independent_thermal(e.field,e.paths,hits,e.metadata['source_RZ_m'],e.thermal_paths.qu,e.reference)
    good=(hits.segment>=0)&np.isfinite(result['thermal']['sigma_m'])
    errors={key:float(np.max(abs(independent[key][good]-result['thermal'][key][good])/independent[key][good]))
            for key in independent}
    row['independent_thermal_max_relative_errors']=errors
    row['audit_passed']=bool(row['audit_passed'] and all(x<=.005 for x in errors.values())
        and row['thermal']['conduction_relative_residual']<=1e-10 and row['thermal']['sheath_relative_residual']<=1e-10)
    return row


def validate_source(folder,levels=LEVELS):
    folder=Path(folder)
    if (folder/'validation.json').exists() or any(folder.glob('*/*_refined.npz')):
        raise FileExistsError(f'existing validation evidence: {folder}')
    verify(ROOT,folder);manifest=load_json(folder/'manifest.json')
    if manifest.get('schema')!='reference-transport-study.v1':raise ValueError('wrong study schema')
    reference=ReferencePoint(**manifest['reference'])
    record=dict(schema='reference-transport-validation.v1',complete=False,levels=[],comparisons={},
                sensitivity={},physical_validation=False)
    for index,(method,count,step,nodes,wall_step) in enumerate(levels):
        start=time.monotonic()
        direct,meta=prepare_evaluator(ROOT,manifest['source'],count=count,step=step,nodes=nodes,
            method=method,bound=manifest['settings']['bound_mm']/1000,boundary_mode='continuous')
        if any(meta[k]!=manifest['metadata'][k] for k in ['shot_sha256','geometry_sha256','surfaces']):
            raise ValueError('validation inputs differ')
        level=dict(interpolation=method,tubes=count,step_m=step,target_nodes=nodes,wall_step_m=wall_step,models={})
        cached=None
        for mode in manifest['models']:
            e=TransportEvaluator(direct,meta,reference,mode=mode,wall_step_m=wall_step,thermal_paths=cached,
                                 quadrature_order=manifest['settings']['quadrature'])
            cached=e.thermal_paths
            base=audit(e,e.baseline)
            rows={'baseline':base}
            for name,c in manifest['candidates'].items():
                if name=='baseline':continue
                result=e.evaluate(c);rows[name]=audit(e,result,base['audit_passed'])
                if index==2:export(folder/mode/(name+'_refined.npz'),e,result)
            if index==2:export(folder/mode/'baseline_refined.npz',e,e.baseline)
            level['models'][mode]=rows
            print(manifest['source'],method,count,mode,'baseline audit',base['audit_passed'],
                  'thermal errors',base.get('independent_thermal_max_relative_errors'),flush=True)
        level['elapsed_seconds']=time.monotonic()-start
        record['levels'].append(level);save(folder/'validation.json',record)
        if index==1:
            variants={
                'half_power':replace(reference,power_W=reference.power_W/2),
                'double_power':replace(reference,power_W=reference.power_W*2),
                'half_density':replace(reference,density_m3=reference.density_m3/2),
                'double_density':replace(reference,density_m3=reference.density_m3*2),
                'half_chi':replace(reference,chi_perp_m2_s=reference.chi_perp_m2_s/2),
                'double_chi':replace(reference,chi_perp_m2_s=reference.chi_perp_m2_s*2)}
            selected={k:v for k,v in manifest['candidates'].items() if k in ['direct_best','thermal_best','frozen_transport_best']}
            for label,ref in variants.items():
                e=TransportEvaluator(direct,meta,ref,mode='thermal',wall_step_m=wall_step,
                                     quadrature_order=manifest['settings']['quadrature'])
                record['sensitivity'][label]={'baseline':summarize(e,e.baseline),
                    'candidates':{k:summarize(e,e.evaluate(c)) for k,c in selected.items()},
                    'scope':'fixed candidate; single-resolution sensitivity, not a new optimum'}
            for anchor in [.0025,.01]:
                e=TransportEvaluator(direct,meta,reference,mode='thermal',wall_step_m=wall_step,
                                     anchor_m=anchor,thermal_paths=cached)
                record['sensitivity'][f'anchor_{anchor:g}m']={'baseline':summarize(e,e.baseline),
                    'candidates':{k:summarize(e,e.evaluate(c)) for k,c in selected.items()},
                    'scope':'local deposition support assumption, not plasma physics calibration'}
            save(folder/'validation.json',record)
    for mode in manifest['models']:
        record['comparisons'][mode]={}
        for candidate in manifest['candidates']:
            original=bool(manifest['models'][mode]['baseline']['baseline_admissible'] and manifest['cross_evaluation'][mode][candidate]['feasible'])
            verdict=acceptance([v['models'][mode][candidate] for v in record['levels']],original)
            if manifest['smoke_test']:verdict['improvement_confirmed']=False
            record['comparisons'][mode][candidate]=verdict
    record['complete']=True;save(folder/'validation.json',record)
    manifest['validation_status']='complete_see_per_candidate_verdicts';save(folder/'manifest.json',manifest)
    return record


def validate_study(output,source='both'):
    results={}
    for name in ['IMP','OMP'] if source=='both' else [source]:
        results[name]=validate_source(Path(output)/name)['comparisons']
    save(Path(output)/'validation_summary.json',results)
    from verification.plot_reference_transport import write_report
    write_report(output)
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',choices=['IMP','OMP','both'],default='both')
    args=parser.parse_args();validate_study(args.output,args.source)
