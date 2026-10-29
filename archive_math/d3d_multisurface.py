"""Explicit DIII-D design components; angle scope independent of design scope."""
import numpy as np
import shapely
from shapely.geometry import LineString, Polygon
from .frozen_field import exponential_source
from .frozen_mapping import Hits, first_hits, map_heat, deform_target, _intersection_points
from .d3d_frozen_case import WIDTHS_MM


class MultiSurfaceEvaluator:
    def __init__(self, field, lcfs, wall, paths, inside, edges, r_sep, z_up,
                 surfaces, target_nodes=161, coefficient_bound=.1, vessel_solids=None):
        self.field, self.lcfs, self.wall = field, np.asarray(lcfs), np.asarray(wall)
        self.paths, self.inside = paths, np.asarray(inside,dtype=bool)
        self.surfaces = [dict(s) for s in surfaces]
        self.bound = float(coefficient_bound)
        self.vessel_solids = vessel_solids
        indices = [s["segment_index"] for s in self.surfaces]
        names = [s["name"] for s in self.surfaces]
        if (not len(indices) or len(set(indices)) != len(indices) or len(set(names)) != len(names)
                or any(type(i) is not int or not 0 <= i < len(wall)-1 for i in indices)
                or self.bound <= 0 or not np.isfinite(self.bound)):
            raise ValueError("unique explicit design surfaces and positive bound required")
        for s in self.surfaces:
            if s.get("angle_floor_deg") not in (None, 1.5):
                raise ValueError("this case permits only an explicit 1.5 deg floor or no floor")
            if "angle_floor_deg" not in s:
                raise ValueError("angle scope must be explicit for each design component")
        self.targets = {s["name"]:np.linspace(wall[s["segment_index"]],wall[s["segment_index"]+1],target_nodes) for s in self.surfaces}
        self.lengths = {n:float(np.linalg.norm(t[-1]-t[0])) for n,t in self.targets.items()}
        self.sources = [exponential_source(field,r_sep,z_up,np.asarray(edges),mm/1000) for mm in WIDTHS_MM]
        self.fixed_indices = np.array([i for i in range(len(wall)-1) if i not in indices])
        self.fixed_geometry = shapely.union_all([LineString(wall[i:i+2]) for i in self.fixed_indices])
        self.fixed_hits = self._fixed_hits()
        self.baseline = self.evaluate(np.zeros(4*len(self.surfaces)))
        if not self.baseline["geometry_valid"]:
            raise ValueError("baseline design components intersect other geometry or LCFS")
        self.baseline_peaks = np.array([float(max(m.q)) for m in self.baseline["maps"]])

    def _fixed_hits(self):
        inter = shapely.intersection(np.asarray(self.paths,dtype=object),self.fixed_geometry)
        ids=np.full(len(self.paths),-1,dtype=int); pts=np.full((len(self.paths),2),np.nan); ds=np.full(len(self.paths),np.inf)
        direction=np.diff(self.wall,axis=0); lengths2=np.sum(direction**2,axis=1)
        for j,(line,g) in enumerate(zip(self.paths,inter)):
            if not self.inside[j]: continue
            candidates=_intersection_points(g)
            if not candidates: continue
            distances=[line.project(shapely.Point(p)) for p in candidates]
            k=int(np.argmin(distances)); point=candidates[k]
            a=np.clip(np.sum((point-self.wall[:-1])*direction,axis=1)/lengths2,0,1)
            residual=np.linalg.norm(self.wall[:-1]+a[:,None]*direction-point,axis=1)
            ids[j]=self.fixed_indices[np.argmin(residual[self.fixed_indices])]
            pts[j]=point; ds[j]=distances[k]
        return Hits(ids,pts,ds)

    def hit_masks(self,result):
        heat=result["maps"][0]
        power=np.max([m.hit_power for m in result["maps"]],axis=0)
        design=np.zeros(len(heat.segment),dtype=bool); angle=design.copy()
        for s in self.surfaces:
            selected=np.isin(heat.segment,result["surface_segments"][s["name"]]) & (power>0)
            design |= selected
            if s["angle_floor_deg"] is not None: angle |= selected
        return design,angle

    def evaluate(self,coefficients):
        coefficients=np.asarray(coefficients,dtype=float)
        if (coefficients.shape != (4*len(self.surfaces),) or not np.all(np.isfinite(coefficients))
                or np.max(np.abs(coefficients)) > self.bound+1e-12):
            raise ValueError("finite coefficient vector must respect the explicit bound")
        targets={s["name"]:deform_target(self.targets[s["name"]],coefficients[4*j:4*j+4]) for j,s in enumerate(self.surfaces)}
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
