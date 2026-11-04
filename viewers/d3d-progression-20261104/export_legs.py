"""Display-only poloidal separatrix legs: fixed psi, first intersection per wall.

No optimization, thermal solve, helical orbit, or L_parallel calculation.
"""
from export_data import HERE,REPO,read,sha,save
import numpy as np
from shapely.geometry import LineString,Point
from shapely.ops import substring
from d3d_optimization.frozen_field import FrozenField
from d3d_optimization.continuous_boundary import extended_leg
from d3d_optimization.frozen_mapping import first_hits

def clip_leg(path,wall,tolerance=1e-6):
    line=path if isinstance(path,LineString) else LineString(path)
    hit=first_hits([line],wall)
    if hit.segment[0]<0:
        return {'resolved':False,'points':[],'strike_RZ_m':None,'reason':'no_first_intersection_with_display_wall'}
    end=hit.points[0];distance=float(hit.distance[0]) if hasattr(hit,'distance') else float(line.project(Point(end)))
    clipped=substring(line,0,distance)
    if not isinstance(clipped,LineString):
        return {'resolved':False,'points':[],'strike_RZ_m':None,'reason':'wall_at_leg_origin'}
    points=np.asarray(clipped.simplify(tolerance,preserve_topology=True).coords)
    if LineString(wall).distance(Point(end))>1e-9:raise ValueError('strike is not on displayed wall')
    return {'resolved':True,'points':points.tolist(),'strike_RZ_m':end.tolist(),
            'wall_segment':int(hit.segment[0]),'poloidal_length_m':distance}

def run():
    review_path=REPO/'baseline_review/baseline_review.json';field_path=REPO/'baseline_review/field.npz';data_path=HERE/'dist/data.json'
    inputs={str(p):sha(p) for p in [review_path,field_path,data_path]}
    data=read(data_path);review=read(review_path)
    for name in ['wall','lcfs']:
        if LineString(data['baseline'][name]).hausdorff_distance(LineString(review[name]))>1e-10:
            raise ValueError('baseline reference geometry mismatch: '+name)
    with np.load(field_path) as z:
        # extended_leg uses psi and its poloidal derivatives only. No B_phi is used.
        field=FrozenField(z['R'],z['Z'],z['psi'],tf_current_ka=0,method='cubic')
    xp=np.asarray(review['xpoint_RZ_m']);psi_sep=float(field.values([xp])[0])
    axis=np.asarray(review['axis_RZ_m']);span=abs(psi_sep-float(field.values([axis])[0]))
    paths={};path_audit={}
    for name,source in [('inner','IMP'),('outer','OMP')]:
        traced=extended_leg(field,xp,psi_sep,source,spacing=.00025)
        points=np.vstack([xp,np.asarray(traced.coords)])
        paths[name]=LineString(points)
        residual=float(np.max(abs(field.values(points)-psi_sep))/span)
        if residual>1e-6:raise ValueError('separatrix psiN drift exceeds display tolerance')
        path_audit[name]={'nodes':len(points),'max_psiN_residual':residual}
    walls={'baseline':data['baseline']['wall']}
    for series in data['series']:
        for row in series['stages']+([series['final']] if series['final'] else []):
            if row['id'] in walls and not np.array_equal(row['wall'],walls[row['id']]):raise ValueError('same identity has differing walls')
            walls[row['id']]=row['wall']
    cases={};unresolved=[]
    for identity,wall in walls.items():
        cases[identity]={name:clip_leg(path,wall) for name,path in paths.items()}
        for name,leg in cases[identity].items():
            if not leg['resolved']:unresolved.append([identity,name]);continue
            assert np.linalg.norm(np.asarray(leg['points'][0])-xp)<1e-12
            assert np.linalg.norm(np.asarray(leg['points'][-1])-leg['strike_RZ_m'])<1e-12
    baseline_errors={name:float(np.linalg.norm(np.asarray(cases['baseline'][name]['strike_RZ_m'])-review['legs'][name]['strike_RZ_m'])) for name in paths}
    if max(baseline_errors.values())>5e-6:raise ValueError('baseline leg strike mismatch')
    for p,h in inputs.items():
        if sha(p)!=h:raise ValueError('protected viewer/research input changed')
    payload={'schema':'d3d-viewer-poloidal-legs.v1','xpoint_RZ_m':xp.tolist(),'cases':cases,
             'data_sha256':inputs[str(data_path)],'scope':'poloidal separatrix from X-point to first hit on displayed wall; axisymmetric display',
             'length_definition':'poloidal arc length only; not upstream parallel connection length',
             'simplification_tolerance_m':1e-6}
    save(HERE/'dist/legs.json',payload)
    save(HERE/'legs-audit.json',{'source_sha256':inputs,'path_checks':path_audit,'baseline_strike_error_m':baseline_errors,
         'cases':len(cases),'unresolved':unresolved,'original_data_unchanged':True,'new_thermal_solves':0,'new_optimizations':0,
         'new_geometric_separatrix_traces':2,'display_leg_file_sha256':sha(HERE/'dist/legs.json')})
    print('Legs:',len(cases),'walls, unresolved:',len(unresolved),'baseline strike error m:',baseline_errors)

if __name__=='__main__':run()
