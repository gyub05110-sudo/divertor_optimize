"""Physical source-domain bounds; a rectangular EFIT grid is not a plasma domain."""
import numpy as np
from shapely.geometry import Polygon,Point,LineString

def physical_source_limit(wall,r_sep,z_up,r_grid_max):
    polygon=Polygon(wall)
    if not polygon.contains(Point(r_sep,z_up)):raise ValueError('Separatrix source lies outside public wall')
    extent=max(float(np.max(np.asarray(wall)[:,0])),r_grid_max)+1
    cross=polygon.intersection(LineString([(r_sep,z_up),(extent,z_up)]))
    parts=[cross] if cross.geom_type=='LineString' else list(cross.geoms)
    segment=next((p for p in parts if p.geom_type=='LineString' and p.distance(Point(r_sep,z_up))<1e-10),None)
    if segment is None:raise ValueError('Cannot identify contiguous outboard source interval')
    wall_r=max(p[0] for p in segment.coords)
    return {'wall_R_m':wall_r,'r_sep_m':r_sep,'z_up_m':z_up,'wall_max_offset_m':wall_r-r_sep,
            'grid_max_offset_m':r_grid_max-r_sep,'physical_max_offset_m':min(wall_r,r_grid_max)-r_sep}

def require_physical_extension(requested_m,physical_max_m):
    if not np.isfinite(requested_m) or requested_m<=0:raise ValueError('Positive source range required')
    if requested_m>physical_max_m+1e-10:raise ValueError('Requested source extends beyond the public wall or EFIT grid')

def scope_reason(psi,psi_start,psi_end):
    t=(np.asarray(psi)-psi_start)/(psi_end-psi_start)
    if t<0:return 'separatrix_side'
    if t>1:return 'beyond_source_end'
    return 'within_source_range'
