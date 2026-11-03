import numpy as np
import pytest
from shapely.geometry import Point, Polygon, LineString

from cpd_optimization.transport import vacuum_boundaries, first_loop_hits, deposit_loops


def test_detached_solid_creates_actual_hole_without_gap_bridge():
    cavity = Polygon([(1,-2),(5,-2),(5,2),(1,2)])
    island = Polygon([(2,-1),(3,-1),(3,0),(2,0)])
    region, loops = vacuum_boundaries(cavity, [island], [4,1])
    assert len(loops) == 2
    assert region.contains(Point(2.5,.1))
    paths = [LineString([(1.5,1),(1.5,-3)]), LineString([(2.5,1),(2.5,-3)])]
    hits, assignment = first_loop_hits(paths, loops)
    np.testing.assert_allclose(hits.points, [[1.5,-2],[2.5,0]])
    assert assignment[0] != assignment[1]
    np.testing.assert_allclose(hits.distance, [3,1])


def test_seed_inside_material_rejected():
    with pytest.raises(ValueError, match='vacuum'):
        vacuum_boundaries(Polygon([(1,0),(4,0),(4,3),(1,3)]),
                          [Polygon([(2,1),(3,1),(3,2),(2,2)])], [2.5,1.5])


class Field:
    def gradient(self, p):
        return np.tile([1.,0.], (len(p),1))


def test_multiple_loops_keep_global_unhit_and_power_without_renormalization():
    loops = [np.array([[1,0],[2,0],[2,1],[1,1],[1,0]],float),
             np.array([[3,0],[4,0],[4,1],[3,1],[3,0]],float)]
    paths = [LineString([(1.5,2),(1.5,.5)]),
             LineString([(3.5,2),(3.5,.5)]),
             LineString([(5,2),(5,.5)])]
    hits, assignment = first_loop_hits(paths, loops)
    result = deposit_loops(Field(), hits, assignment, loops,
                           np.array([.2,.3,.5]), np.array([.1,.1,np.nan]), 1., .05)
    assert result['ledger']['prescribed'] == pytest.approx(1.)
    assert result['ledger']['deposited'] == pytest.approx(.5)
    assert result['ledger']['unhit'] == pytest.approx(.5)
    assert result['ledger']['balance_error'] < 1e-12
    assert [d['ledger']['deposited'] for d in result['depositions']] == pytest.approx([.2,.3])
    assert result['peak_q'] == max(np.max(d['q']) for d in result['depositions'])
    assert result['area_90_m2'] > 0


def test_nonfinite_total_power_rejected_even_when_all_donors_unhit():
    loops=[np.array([[1,0],[2,0],[2,1],[1,1],[1,0]],float)]
    hits, assignment=first_loop_hits([LineString([(3,0),(3,1)])]*2,loops)
    with np.errstate(over='ignore'), pytest.raises(ValueError):
        deposit_loops(Field(),hits,assignment,loops,np.array([1e308,1e308]),
                      np.array([np.nan,np.nan]),1.,.05)


def test_loop_order_and_orientation_preserve_earliest_hit():
    loops=[np.array([[1,0],[2,0],[2,1],[1,1],[1,0]],float),
           np.array([[3,0],[4,0],[4,1],[3,1],[3,0]],float)]
    paths=[LineString([(5,.5),(.5,.5)])]
    first,_=first_loop_hits(paths,loops)
    second,_=first_loop_hits(paths,[a[::-1] for a in loops[::-1]])
    np.testing.assert_array_equal(first.points,second.points)
    np.testing.assert_array_equal(first.distance,second.distance)
