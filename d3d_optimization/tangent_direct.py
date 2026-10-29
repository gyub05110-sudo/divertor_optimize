"""Frozen direct-map and witness adapters: only curve construction differs.

Copied from multisurface.MultiSurfaceEvaluator.evaluate and
continuous_boundary.BoundaryGuard.witness on 2026-10-29 so archived files remain
byte-identical. No runtime monkeypatching or global shape substitution.
"""
import numpy as np
import shapely
from shapely.geometry import LineString,Polygon
from scipy.optimize import brentq
from .frozen_mapping import first_hits,map_heat,Hits
from .shape import ReferencePolyline
from .multisurface import MultiSurfaceEvaluator
from .continuous_boundary import BoundaryGuard
from .tangent_geometry import Curve

class TangentMultiSurface(MultiSurfaceEvaluator):
    def curve(self,name,c):
        a,b=self.targets[name][[0,-1]]
        direction=None
        if name=='floor':
            surface=next(s for s in self.surfaces if s['name']==name)
            i=surface['segment_index']
            direction=self.wall[i+2]-self.wall[i+1]
        return Curve(a,b,c,direction)

    def evaluate(self,coefficients):
        coefficients=np.asarray(coefficients,dtype=float)
        if (coefficients.shape != (4*len(self.surfaces),) or not np.all(np.isfinite(coefficients))
                or np.max(np.abs(coefficients)) > self.bound+1e-12):
            raise ValueError("finite coefficient vector must respect the explicit bound")
        targets={s["name"]:self.curve(s["name"],coefficients[4*j:4*j+4]).point(np.linspace(0,1,len(self.targets[s["name"]]))) for j,s in enumerate(self.surfaces)}
        invalid={"geometry_valid":False,"feasible":False,"penalty":1000.,"targets":targets}
        for name,t in targets.items():
            if (np.sum(np.linalg.norm(np.diff(t,axis=0),axis=1)) > 1.2*self.lengths[name]
                    or not LineString(t).is_simple or LineString(t[1:-1]).intersects(self.fixed_geometry)
                    or LineString(t).intersects(Polygon(self.lcfs))):
                return invalid
            if self.vessel_solids is not None and LineString(t).intersects(self.vessel_solids):
                return invalid
        values=list(targets.values())
        for i in range(len(values)):
            for j in range(i):
                if LineString(values[i][1:-1]).intersects(LineString(values[j][1:-1])): return invalid
        by_index={s["segment_index"]:s["name"] for s in self.surfaces}
        vertices=[self.wall[0]]; original=[]; surface_segments={}
        for i in range(len(self.wall)-1):
            if i in by_index:
                name=by_index[i]; start=len(original); t=targets[name]
                vertices.extend(t[1:]); original.extend([i]*(len(t)-1))
                surface_segments[name]=np.arange(start,len(original))
            else:
                vertices.append(self.wall[i+1]); original.append(i)
        candidate=np.asarray(vertices); original=np.asarray(original,dtype=int)
        orig_to_candidate={int(orig):int(k) for k,orig in enumerate(original) if orig in self.fixed_indices}
        segment=np.array([orig_to_candidate.get(int(i),-1) for i in self.fixed_hits.segment])
        points=self.fixed_hits.points.copy(); distances=self.fixed_hits.distance.copy()
        for s in self.surfaces:
            name=s["name"]; hit=first_hits(self.paths,targets[name])
            choose=self.inside & (hit.distance < distances)
            segment[choose]=surface_segments[name][0]+hit.segment[choose]
            points[choose]=hit.points[choose]; distances[choose]=hit.distance[choose]
        maps=[map_heat(self.field,Hits(segment,points,distances),src.weights,src.h,candidate) for src in self.sources]
        active=np.max([m.hit_power for m in maps],axis=0)
        angles={}; angle_penalty=0.; angle_ok=True
        for s in self.surfaces:
            name=s["name"]; floor=s["angle_floor_deg"]
            status=maps[0].angle_application(floor if floor is not None else 0.,surface_segments[name],active)
            status["angle_floor_deg"]=floor
            if floor is None:
                status.update(status="diagnostic_only_no_angle_floor",angle_constraint_passed=None,applied_hit_count=0)
            elif status["status"] == "target_angle_violation":
                angle_ok=False; angle_penalty+=100*(floor-status["minimum_angle_deg"])
            elif status["status"] == "invalid_target_angle":
                angle_ok=False; angle_penalty+=1000
            # An unheated shelf has no angle violation. The overall design must
            # still intercept power; transferring heat to other walls is counted.
            angles[name]=status
        power={n:[float(np.sum(m.segment_power[ids])) for m in maps] for n,ids in surface_segments.items()}
        missing=max(m.unmapped_power for m in maps); balance=max(m.balance_error for m in maps)
        illuminated=any(any(v>0 for v in ps) for ps in power.values())
        feasible=bool(angle_ok and illuminated and missing<=1e-3 and balance<=1e-4)
        penalty=angle_penalty+1000*max(0,missing-1e-3)+(0 if illuminated else 1000)
        return {"geometry_valid":True,"feasible":feasible,"penalty":penalty,"targets":targets,
                "wall":candidate,"original_segment":original,"surface_segments":surface_segments,
                "maps":maps,"angles":angles,"surface_power":power,"missing":missing,
                "angle_penalty":angle_penalty,"power_balance_error":balance}


class TangentGuard(BoundaryGuard):
    def witness(self, result, coefficients):
        hit=first_hits([self.path],result['wall'])
        k=int(hit.segment[0])
        if k<0: return dict(resolved=False,reason='no_separatrix_first_hit')
        original=int(result['original_segment'][k])
        surface=next(((j,s) for j,s in enumerate(self.e.surfaces) if s['segment_index']==original),None)
        if surface:
            j,surface=surface
            target=self.e.targets[surface['name']]
            ref=ReferencePolyline([target[0],target[-1]])
            spline=self.e.curve(surface["name"],coefficients[4*j:4*j+4])
            point=lambda s: spline.point(s)
            tangent=lambda s: spline.tangent(s)
            ids=result['surface_segments'][surface['name']]
            index=k-int(ids[0]); n=len(ids)
            lo=max(0,(index-2)/n); hi=min(1,(index+3)/n)
        else:
            a,b=result['wall'][k:k+2]
            point=lambda s: a+s*(b-a)
            tangent=lambda s: b-a
            lo,hi=0.,1.
        def residual(s): return float(self.e.field.values([point(s)])[0]-self.psi_sep)
        grid=np.linspace(lo,hi,13)
        vals=[residual(s) for s in grid]
        roots=[float(s) for s,v in zip(grid,vals) if abs(v)<1e-14]
        roots += [brentq(residual,l,h,xtol=1e-14) for l,h,u,v in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]) if u*v<0]
        if not roots:return dict(resolved=False,reason='first_hit_flux_root_not_bracketed')
        s=min(roots,key=lambda s:np.linalg.norm(point(s)-hit.points[0]))
        p=point(s); t=tangent(s); normal=np.array([-t[1],t[0]])/np.linalg.norm(t)
        delta=float(np.linalg.norm(p-hit.points[0]))
        # Local root refinement must remain near the observed first intersection.
        if delta>2e-3 or abs(residual(s))>1e-11:
            return dict(resolved=False,reason='first_hit_root_refinement_moved_too_far')
        magnetic=self.e.field.magnetic([p])[0]
        q=self.h0*abs(np.dot(magnetic[:2],normal))
        angle=float(np.degrees(np.arcsin(np.clip(abs(np.dot(magnetic[:2],normal))/np.linalg.norm(magnetic),0,1))))
        floor=surface.get('angle_floor_deg') if surface else None
        return dict(resolved=True,point=p.tolist(),q=q.tolist(),angle_deg=angle,
                    angle_passed=bool(floor is None or angle>=floor),surface=surface['name'] if surface else 'fixed',
                    split_segment=original,psi_residual=residual(s),refinement_distance_m=delta,
                    normal_RZ=normal.tolist(),zero_power_limit=True)
