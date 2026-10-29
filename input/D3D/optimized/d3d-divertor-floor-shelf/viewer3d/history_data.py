"""Scientific display helpers. Curves are derived from the frozen EFIT, not sketches."""
import numpy as np
from shapely.geometry import Polygon,Point
from shapely.ops import substring

def select_history_records(records,safe):
    baseline=next(r for r in records if r['stage']=='baseline')
    result=[baseline];best=float(baseline['worst_all_wall_peak_ratio']) if baseline['feasible'] else float('inf')
    for row in records:
        if row['case']==baseline['case']:continue
        if row['feasible'] and safe.get(row['case'],False) and row['worst_all_wall_peak_ratio']<best-1e-12:
            result.append(row);best=float(row['worst_all_wall_peak_ratio'])
    return result

def normalize_psi(psi,psi_axis,psi_sep):
    span=float(psi_sep-psi_axis)
    if not np.isfinite(span) or abs(span)<1e-12:raise ValueError('distinct finite axis and separatrix psi required')
    return (np.asarray(psi)-psi_axis)/span

def truncate_path(path,hit,stride=8):
    point=Point(hit);distance=path.project(point)
    if path.distance(point)>1e-7:raise ValueError('first hit must lie on its SOL trace')
    shortened=substring(path,0,distance)
    coords=np.asarray(shortened.coords)
    sampled=coords[::stride]
    if not np.allclose(sampled[-1],coords[-1],atol=1e-12):sampled=np.vstack([sampled,coords[-1]])
    return sampled.tolist()

def core_surfaces(field,axis,psi_sep):
    import contourpy
    psi_axis=float(field.values([axis])[0])
    rr=np.linspace(field.r[0],field.r[-1],240);zz=np.linspace(field.z[0],field.z[-1],320)
    rg,zg=np.meshgrid(rr,zz)
    normalized=normalize_psi(field.values(np.column_stack([rg.ravel(),zg.ravel()])).reshape(rg.shape),psi_axis,psi_sep)
    contour=contourpy.contour_generator(x=rr,y=zz,z=normalized)
    result=[]
    for level in [.2,.4,.6,.8,.95]:
        closed=[line for line in contour.lines(level) if len(line)>8 and np.linalg.norm(line[-1]-line[0])<1e-8
                and Polygon(line).contains(Point(axis))]
        if not closed:raise ValueError(f'no closed axis-enclosing flux contour at psi_N={level}')
        curve=max(closed,key=lambda line:Polygon(line).area)
        residual=float(np.max(np.abs(normalize_psi(field.values(curve),psi_axis,psi_sep)-level)))
        if residual>.01:raise ValueError('contour interpolation residual exceeds visualization tolerance')
        result.append({'psi_normalized':level,'points':curve.tolist(),'max_psiN_residual':residual})
    return result,{'psi_axis':psi_axis,'psi_separatrix':float(psi_sep),'axis_RZ_m':axis.tolist(),
                   'definition':'(psi - psi_axis)/(psi_Xpoint - psi_axis)','method':'cubic EFIT spline; closed axis-enclosing contours'}
