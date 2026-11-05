"""Read dimensional heat evidence at the actual source power of its physics solve."""
import numpy as np


def read_physical_profile(path, operating_power_W=None):
    with np.load(path,allow_pickle=False) as z:
        recorded=float(z['assigned_power_W']) if 'assigned_power_W' in z else None
        if recorded is not None and operating_power_W is not None and recorded!=operating_power_W:
            raise ValueError('operating power mismatch with stored physics run')
        power=recorded if recorded is not None else operating_power_W
        if power is None or not np.isfinite(power) or power<=0:
            raise ValueError('provide the recorded operating power in W for a legacy export')
        wall=z['receiving_wall'];area=np.diff(z['area_edges']);fraction=z['cell_power']
        points=z['points'];q=z['q_MW_m2']
        if (wall.shape!=(len(area)+1,2) or points.shape!=(q.size,2) or q.ndim!=1 or fraction.shape!=area.shape
                or not all(np.isfinite(v).all() for v in (wall,area,fraction,points,q))
                or np.any(area<=0) or np.any(fraction<0) or np.any(q<0) or np.any(wall[:,0]<=0)
                or np.linalg.norm(wall[-1]-wall[0])>1e-10):
            raise ValueError('invalid physical wall or heat profile')
        physical_area=np.pi*(wall[:-1,0]+wall[1:,0])*np.linalg.norm(np.diff(wall,axis=0),axis=1)
        if not np.allclose(area,physical_area,rtol=2e-8,atol=1e-13):
            raise ValueError('area units or physical geometry mismatch')
        cell_MW=fraction*(power/1e6)
        return dict(assigned_power_W=power,wall=wall,area_m2=area,cell_power_MW=cell_MW,
                    cell_q_MW_m2=cell_MW/area,points=points,continuous_q_MW_m2=q)

