import json
from pathlib import Path
import numpy as np
import pytest
from shapely.geometry import LineString,Point
from d3d_optimization.frozen_field import FrozenField,exponential_source
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.independent_sources import midplane_source,split_design_wall,prepare_evaluator

ROOT=Path(__file__).resolve().parents[1]
def test_omp_source_matches_preserved_kernel():
    rr=np.linspace(1,3,33);zz=np.linspace(-2,2,33);r,z=np.meshgrid(rr,zz)
    f=FrozenField(rr,zz,r*r+z*z,100.)
    edges=np.linspace(0,.02,17)
    old=exponential_source(f,2.,0.,edges,.003)
    new=midplane_source(f,2.,0.,edges,.003,'OMP')
    np.testing.assert_array_equal(old.weights,new.weights);np.testing.assert_array_equal(old.h,new.h)

def test_imp_uses_actual_smaller_R_and_normalizes_independently():
    class Field:
        def __init__(self):self.visited=[]
        def magnetic(self,p):
            self.visited.extend(np.array(p)[:,0]);return np.tile([0,1,2.],(len(p),1))
        def gradient(self,p):return np.c_[np.array(p)[:,0]**2,np.zeros(len(p))]
    f=Field();edges=np.linspace(0,.02,17)
    s=midplane_source(f,1.1,0.,edges,.003,'IMP')
    assert max(f.visited)<1.1 and min(f.visited)>1.08
    assert abs(s.weights.sum()-1)<1e-14
    other=midplane_source(Field(),2.25,0.,edges,.003,'OMP')
    assert not np.allclose(s.weights,other.weights,rtol=1e-8,atol=1e-12)

def test_split_wall_keeps_all_excluded_material_and_exact_cutoff():
    d=json.loads((ROOT/'baseline_review/baseline_review.json').read_text())
    w,s,original=split_design_wall(d)
    assert LineString(w).hausdorff_distance(LineString(d['wall']))<1e-12
    assert len(s)==6 and {r['original_segment'] for r in s}=={58,59,60,68,69,70}
    for surface in s:
        if surface['original_segment'] in [69,70]:
            assert max(np.array(surface['endpoints_RZ_m'])[:,1])<=-1.15
    fixed=set(range(len(w)-1))-{r['segment_index'] for r in s}
    assert {53,54,55}.issubset({original[i] for i in fixed})

@pytest.mark.parametrize('source',['IMP','OMP'])
def test_actual_source_cached_hits_and_fixed_wall_survive_deformation(source):
    e,meta=prepare_evaluator(ROOT,source,count=48,step=.002,nodes=41)
    c=np.full(e.dimension,.0002);result=e.evaluate(c)
    assert result['geometry_valid']
    full=first_hits(e.paths,result['wall'])
    np.testing.assert_array_equal(full.segment,result['maps'][0].segment)
    np.testing.assert_allclose(full.points,result['maps'][0].points,atol=1e-12)
    for m in result['maps']:assert abs(m.mapped_power+m.unmapped_power-1)<1e-14
    for i in e.fixed_indices:
        assert LineString(result['wall']).distance(Point(e.wall[i]))<1e-12
    assert meta['source']==source and meta['source_power']==1
    assert all(abs(s.weights.sum()-1)<1e-14 for s in e.sources)
