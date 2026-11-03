"""CPD fixed-equilibrium adapter; keeps existing D3D physics modules unchanged."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import hashlib
import json
import numpy as np
import shapely
from scipy.optimize import brentq
from shapely.geometry import LineString, Point, Polygon

from .geometry import load_geometry
from .equilibrium import FrozenField
from .transport import vacuum_boundaries, first_loop_hits, deposit_loops
from d3d_optimization.independent_sources import midplane_source
from d3d_optimization.frozen_mapping import trace_paths
from d3d_optimization.reference_transport import ReferencePoint, ThermalPaths, thermal_summary


def deadline_remaining(start, seconds, *, now=None):
    return max(0., float(seconds)-((perf_counter() if now is None else now)-start))


def midplane_separatrix(field, axis, psi_sep, source):
    if source not in ('IMP', 'OMP'):
        raise ValueError('source must be IMP or OMP')
    stop = field.r[0]+1e-8 if source == 'IMP' else field.r[-1]-1e-8
    radii = np.linspace(axis[0], stop, 2001)
    values = field.values(np.c_[radii, np.full(len(radii), axis[1])])-psi_sep
    for a,b,fa,fb in zip(radii[:-1],radii[1:],values[:-1],values[1:]):
        if fa*fb < 0:
            return float(brentq(lambda r: float(field.values([[r,axis[1]]])[0]-psi_sep),
                                min(a,b),max(a,b),xtol=1e-13))
    raise ValueError(f'No {source} separatrix crossing from selected magnetic axis')


class CPDEvaluator:
    def __init__(self, nodes, equilibrium, source, *, count=160, path_step=.002,
                 wall_step=.003, arc_step=.003, max_path=30., quadrature=64,
                 cache_dir=None, allow_diagnostic=False):
        started=perf_counter()
        self.nodes=Path(nodes);self.equilibrium=Path(equilibrium)
        report=json.loads(self.equilibrium.with_suffix('.json').read_text(encoding='utf-8'))
        validated=report.get('status') == 'validated_for_frozen_transport'
        if not validated and not allow_diagnostic:
            raise ValueError('Validated frozen equilibrium required before campaign')
        if count<8 or min(path_step,wall_step,arc_step,max_path)<=0:
            raise ValueError('Positive numerical settings and >=8 donors required')
        self.field=FrozenField(equilibrium)
        if self.field.sha256 != report['frozen_field_sha256']:
            raise ValueError('Frozen equilibrium hash mismatch')
        self.geometry=load_geometry(nodes,max_chord_m=arc_step)
        self.geometry_hash=hashlib.sha256(self.nodes.read_bytes()).hexdigest()
        if self.geometry_hash != report['source_sha256']:
            raise ValueError('Equilibrium geometry differs from baseline input')
        with np.load(equilibrium,allow_pickle=False) as data:
            self.axis=data['axis'].copy();self.lcfs=data['lcfs'].copy()
        self.cavity=Polygon(self.geometry.contours[('VV','inner')].points_m)
        self.vacuum,self.loops=vacuum_boundaries(self.cavity,
            [self.geometry.solids[k] for k in ('Blanket','Shield','Divertor')],self.axis)
        div=self.geometry.solids['Divertor']
        matches=[k for k,loop in enumerate(self.loops)
                 if LineString(loop).hausdorff_distance(div.boundary)<1e-8]
        if len(matches)!=1:
            raise ValueError('Detached divertor must correspond to one entire vacuum loop')
        self.divertor_loop=matches[0]
        faces=self.geometry.contours[('Divertor','outer')].faces
        tree=shapely.STRtree([face.line for face in faces])
        wall=self.loops[self.divertor_loop]
        ids=tree.nearest(shapely.points((wall[:-1]+wall[1:])/2))
        self.base_face_ids=np.array([faces[k].start_node_id for k in ids],int)
        self.source_name=source;self.reference=ReferencePoint();self.wall_step=wall_step
        radius=midplane_separatrix(self.field,self.axis,self.field.psi_lcfs,source)
        self.source_rz=[radius,float(self.axis[1])]
        edges=.1*np.linspace(0,1,count+1)**2
        self.source=midplane_source(self.field,*self.source_rz,edges,self.reference.lambda_q_m,source)
        sign=-1 if source=='IMP' else 1
        starts=np.c_[radius+sign*self.source.x,np.full(count,self.axis[1])]
        if not all(self.vacuum.contains(Point(p)) for p in starts):
            raise ValueError('Fixed 100 mm source support leaves connected vacuum')
        self.metadata=dict(source=source,source_RZ_m=self.source_rz,source_support_m=.1,
            count=count,path_step_m=path_step,wall_step_m=wall_step,arc_step_m=arc_step,
            max_path_m=max_path,quadrature=quadrature,equilibrium_sha256=self.field.sha256,
            geometry_sha256=self.geometry_hash,reference=asdict(self.reference),
            independent_source=True,diagnostic_equilibrium=not validated)
        cache_key=hashlib.sha256(json.dumps({k:self.metadata[k] for k in
            ('source','source_RZ_m','count','path_step_m','max_path_m','equilibrium_sha256')},
            sort_keys=True).encode()).hexdigest()
        cache=None if cache_dir is None else Path(cache_dir)/f'paths-{cache_key}.npz'
        reused=bool(cache is not None and cache.exists())
        if reused:
            with np.load(cache,allow_pickle=False) as data:
                xy=data['points'];offsets=data['offsets'];paths=[xy[a:b] for a,b in zip(offsets[:-1],offsets[1:])]
        else:
            paths=trace_paths(self.field,starts,step=path_step,max_length=max_path,
                              extend_below=self.field.z[0]+2*path_step)
            if cache is not None:
                cache.parent.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(cache,points=np.concatenate(paths),offsets=np.r_[0,np.cumsum([len(p) for p in paths])])
        self.paths=[LineString(p) for p in paths]
        self.thermal_paths=ThermalPaths(self.field,self.paths,self.source_rz,self.source,
                                       self.reference,quadrature_order=quadrature)
        self.metadata['path_cache_reused']=reused
        self.metadata['preparation_seconds']=perf_counter()-started
        self._deformation=None

    @property
    def deformation(self):
        if self._deformation is None:
            from .deformation import DivertorDeformation
            self._deformation=DivertorDeformation(self.geometry,lcfs=Polygon(self.lcfs),
                                                 max_chord_m=self.metadata['arc_step_m'])
        return self._deformation

    def evaluate(self,coefficients=None):
        started=perf_counter();loops=list(self.loops)
        geometry_metrics=None
        face_ids=self.base_face_ids
        if coefficients is not None:
            geometry_metrics=self.deformation.evaluate(coefficients)
            if not geometry_metrics['feasible']:
                return dict(feasible=False,errors=geometry_metrics['errors'],elapsed_seconds=perf_counter()-started)
            loops[self.divertor_loop]=geometry_metrics['wall']
            face_ids=geometry_metrics['segment_face_ids']
        hits,assignment=first_loop_hits(self.paths,loops)
        thermal=self.thermal_paths.evaluate(hits)
        regions=[np.zeros(len(loop)-1,int) for loop in loops]
        regions[self.divertor_loop]=np.where((face_ids>=6)&(face_ids<=17),1,2)
        dep=deposit_loops(self.field,hits,assignment,loops,self.source.weights,
                          thermal['sigma_m'],self.thermal_paths.grad_sep,self.wall_step,regions)
        result=dict(feasible=dep['ledger']['unhit']<=.001 and dep['peak_q']>0,
            errors=[] if dep['ledger']['unhit']<=.001 else ['unhit_power'],
            loops=loops,hits=hits,assignment=assignment,thermal=thermal,deposition=dep,
            geometry=geometry_metrics,divertor_face_ids=face_ids,
            coefficients=None if coefficients is None else np.asarray(coefficients,float).copy(),
            elapsed_seconds=perf_counter()-started)
        return result

    def describe(self,result):
        out={k:result[k] for k in ('feasible','errors','elapsed_seconds')}
        if 'deposition' not in result:return out
        dep=result['deposition'];w=self.source.weights
        peak_dep=max(dep['depositions'],key=lambda d:float(np.max(d['q'])),default=None)
        out.update(peak_MW_m2=dep['peak_q']*self.reference.power_W/1e6,
            peak_RZ_m=None if peak_dep is None else peak_dep['points'][np.argmax(peak_dep['q'])].tolist(),
            ledger=dep['ledger'],area_90_m2=dep['area_90_m2'],effective_area_m2=dep['effective_area_m2'],
            divertor_first_hit_power_fraction=float(w[result['assignment']==self.divertor_loop].sum()),
            thermal=thermal_summary(result['thermal'],w,self.reference),
            per_loop_power={str(d['loop_index']):d['ledger']['deposited'] for d in dep['depositions']},
            diagnostic_equilibrium=self.metadata['diagnostic_equilibrium'])
        div_dep=next((d for d in dep['depositions'] if d['loop_index']==self.divertor_loop),None)
        region_power=np.zeros(3) if div_dep is None else div_dep['donor_region_power'].sum(axis=0)
        out['divertor_front_power_fraction']=float(region_power[1])
        out['divertor_rear_and_top_cap_power_fraction']=float(region_power[2])
        if result['geometry'] is not None:
            out['length_ratios']=result['geometry']['length_ratios']
            out['min_curvature_radius_m']=result['geometry']['min_curvature_radius_m']
        return out

    def export(self,path,result):
        """Retain receiving curves, flux and donor integrals for audit/replay."""
        arrays=dict(hit_points=result['hits'].points,hit_loop=result['assignment'],
                    hit_segment=result['hits'].segment,source_weights=self.source.weights,
                    source_edges=self.source.edges,source_rz=self.source_rz,
                    source_power_W=self.reference.power_W,
                    q_units='m^-2 per unit assigned power; multiply by source_power_W for W/m^2',
                    equilibrium_sha256=self.field.sha256,geometry_sha256=self.geometry_hash)
        arrays['divertor_loop_index']=self.divertor_loop
        arrays['divertor_segment_start_node_ids']=result['divertor_face_ids']
        arrays['coefficients']=np.array([]) if result['coefficients'] is None else result['coefficients']
        arrays.update({'thermal_'+key:value for key,value in result['thermal'].items()})
        for k,loop in enumerate(result['loops']):arrays[f'loop_{k}']=loop
        for d in result['deposition']['depositions']:
            k=d['loop_index']
            for name in ('points','q','cell_power','areas','wall','cell_parent','sample_parent'):
                arrays[f'loop_{k}_{name}']=d[name]
        np.savez_compressed(path,**arrays)
