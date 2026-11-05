from dataclasses import asdict
import numpy as np
import pytest
from d3d_optimization.reference_transport import ReferencePoint
from verification.nf_search_benchmark import triangle, shortlist
from run_global_campaign import campaign_reference

def test_triangle_includes_each_anchor_once_and_rejects_collinear():
    p=np.array([[0.,0.],[1.,0.],[0.,1.]])
    rows=triangle(p,6)
    assert len(rows)==28 and len({tuple(row) for row in rows})==28
    for anchor in p:
        assert np.any(np.all(rows==anchor,axis=1))
    with pytest.raises(ValueError):triangle(np.zeros((3,2)),6)

def test_shortlist_rejects_infeasible_and_handles_ties_and_empty():
    rows=[{'i':i,'eligible':ok,'peak':v} for i,ok,v in [(4,False,.01),(2,True,.8),(1,True,.8),(9,True,1.2)]]
    assert shortlist(rows,2)==[1,2] and shortlist([],3)==[]
    with pytest.raises(ValueError):shortlist([{'i':1,'eligible':True,'peak':float('nan')}],2)

def test_campaign_reference_converts_units_once_and_preserves_defaults():
    assert campaign_reference({})==asdict(ReferencePoint())
    x=campaign_reference({'power_mw':2.,'density':3e19,'chi':.6,'lambda_q_mm':1.5,'gamma':8.,'kappa0':1600.})
    assert x==asdict(ReferencePoint(2e6,3e19,.6,.0015,8.,1600.))
    with pytest.raises(ValueError):campaign_reference({'chi':0.})


def test_launch_cannot_override_frozen_reference(tmp_path):
    import subprocess
    import sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    result=subprocess.run([sys.executable,str(root/'run_global_campaign.py'),'--mode','launch',
                           '--output',str(tmp_path/'absent'),'--chi','0.6'],capture_output=True,text=True)
    assert result.returncode==2 and 'frozen at creation' in result.stderr
    assert not (tmp_path/'absent').exists()

