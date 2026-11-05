import numpy as np
import pytest
from verification.preimpact_routing import Wall, ring_grid, complete_bounds, route_bins, split_routes

class LinearFlux:
    def values(self,p):return np.asarray(p)[:,0]
    def gradient(self,p):return np.tile([1.,0.],(len(p),1))

def test_same_label_on_front_and_back_deposits_once():
    wall=Wall(np.array([[1.,0.],[2.,.5],[2.1,1.],[1.1,1.],[1.,0.]]),LinearFlux())
    grid=ring_grid(wall,.025,0)
    arcs=wall.arc_id(np.array([.2,.3]))
    assert arcs[0]==arcs[1]
    mass=lambda lo,hi: np.asarray(hi)-np.asarray(lo)
    p,cap=route_bins(wall,grid,[(1.2,1.4,int(arcs[0]))],[(1.4,1.5)],mass)
    assert p.sum()==pytest.approx(.2)
    assert p[grid['s_mid']>wall.s[1]].sum()==pytest.approx(0)
    assert cap.sum()>=.1
    assert abs(sum(grid['area'])-wall.area[-1])<1e-11

def test_bounds_include_any_allowed_partial_allocation():
    p=np.array([.1,.6]);a=np.array([.07,.1]);u=.3;cap=np.array([.09,.3])
    b=complete_bounds(p,a,u,cap)
    for i in np.linspace(0,.09,11):
        for j in np.linspace(0,.3-i,11):
            t=p+[i,j]
            peak=max(t/a);A=t.sum()**2/np.sum(t*t/a)
            assert b['peak'][0]-1e-12<=peak<=b['peak'][1]+1e-12
            assert b['A2'][0]-1e-12<=A<=b['A2'][1]+1e-12

def test_disagreement_becomes_unresolved_and_is_never_an_interpolated_hit():
    def label(x): return np.where(x<.31,0,1)
    good,bad=split_routes(np.array([0.,1.]),label,depth=5)
    assert any(lo<=.31<=hi for lo,hi in bad)
    assert sum(h-l for l,h,_ in good)+sum(h-l for l,h in bad)==pytest.approx(1)

