"""Strict, conservative conversion of evaluated cell powers to display bins."""
import numpy as np
from shapely.geometry import LineString

def make_profile(z, expected_wall, expected_peak, power_W):
    wall=np.asarray(z['wall']);mesh=np.asarray(z['receiving_wall'])
    area=np.diff(z['area_edges']);power=np.asarray(z['cell_power']);parent=z['cell_parent']
    q=np.asarray(z['q_MW_m2']);points=np.asarray(z['points'])
    if not all(np.all(np.isfinite(a)) for a in [wall,mesh,area,power,q,points]):
        raise ValueError('nonfinite profile arrays')
    if len(mesh)!=len(power)+1 or len(area)!=len(power) or np.any(area<=0) or np.any(power<0) or np.any(q<0):
        raise ValueError('invalid receiving mesh or negative power')
    if LineString(expected_wall).hausdorff_distance(LineString(wall))>2e-5:
        raise ValueError('candidate/profile wall mismatch')
    if abs(float(power.sum())-1)>2e-12:raise ValueError('deposited power not conserved')
    k=int(np.argmax(q));peak=float(q[k]);relative=abs(peak/expected_peak-1)
    if relative>1e-7:raise ValueError(f'recorded peak mismatch: {peak} vs {expected_peak}')
    segments=[];areas=[];i=0
    while i<len(power):
        j=i+1
        while j<min(i+4,len(power)) and parent[j]==parent[i]:j+=1
        a,b=mesh[i],mesh[j];area_sum=float(area[i:j].sum())
        segments.append([*map(float,a),*map(float,b),float(power[i:j].sum()*power_W/1e6/area_sum)])
        areas.append(area_sum);i=j
    reconstructed=sum(s[4]*a*1e6/power_W for s,a in zip(segments,areas))
    if abs(reconstructed-float(power.sum()))>2e-12:raise ValueError('display-bin power mismatch')
    return dict(segments=segments,areas_m2=areas,peak=peak,peak_point=points[k].tolist(),
        power_fraction=float(power.sum()),display_power_fraction=reconstructed,recorded_peak_relative_error=relative,
        note='Area-weighted average of up to four consecutive cells within each parent segment; full wall coverage. Peak is separately refined within the evaluated grid.')
