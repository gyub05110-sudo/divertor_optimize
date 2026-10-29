"""Regional coefficient restriction without modifying the frozen evaluators."""
import numpy as np

REGIONS={'IMP':('imp_sloped','imp_vertical'),'OMP':('omp_entry','shelf','omp_connector','floor')}


class RegionalEvaluator:
    def __init__(self,full,source):
        self.full=full;self.active_names=REGIONS[source]
        self.indices=np.array([4*j+k for j,s in enumerate(full.surfaces) if s['name'] in self.active_names for k in range(4)])
        if len(self.indices)!=4*len(self.active_names):raise ValueError('missing named design surface')
        self.dimension=len(self.indices)
        self.lengths={n:full.lengths[n] for n in self.active_names}
        self.metadata=dict(full.metadata,active_surfaces=list(self.active_names),active_dimension=self.dimension)

    def __getattr__(self,name):return getattr(self.full,name)

    def expand(self,c):
        c=np.asarray(c,float)
        if c.shape!=(self.dimension,) or not np.all(np.isfinite(c)):raise ValueError('finite regional coefficients required')
        result=np.zeros(self.full.dimension);result[self.indices]=c;return result

    def project(self,c):
        c=np.asarray(c,float)
        if c.shape!=(self.full.dimension,):raise ValueError('full coefficient vector required')
        return c[self.indices].copy()

    def evaluate(self,c):return self.full.evaluate(self.expand(c))

    def describe(self,r):
        row=self.full.describe(r)
        if r.get('geometry_valid'):
            fixed=[s['name'] for s in self.full.surfaces if s['name'] not in self.active_names]
            row['opposite_region_displacement_m']=max((float(np.max(np.linalg.norm(
                r['targets'][n]-self.full.baseline['targets'][n],axis=1))) for n in fixed),default=0.)
        row['active_surfaces']=list(self.active_names)
        return row
