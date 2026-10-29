import importlib.util
import numpy as np
from shapely.geometry import LineString
from d3d_optimization.broadening import PathLengths
from d3d_optimization.frozen_mapping import first_hits
from d3d_optimization.transport_deposition import DepositionControls, deposit_sigma
from verification.test_broadening import LinearField


def test_diagnostic_reconciles_frozen_kernel_without_changing_it():
    assert importlib.util.find_spec('verification.diagnose_boundary_transport') is not None
    from verification.diagnose_boundary_transport import legacy_routes
    field=LinearField();wall=np.array([[1.,0.],[3.,0.],[3.,2.],[1.,2.],[1.,0.]])
    paths=[LineString([[2.,1.],[2.,-.5]])];lengths=PathLengths(field,paths)
    hits=first_hits(paths,wall);cfg=DepositionControls(.05,.005,.001)
    dep=deposit_sigma(field,paths,lengths,hits,wall,[1.],[1.,1.],1,None,[],cfg,[.8])
    row=legacy_routes(paths,hits,dep,[1.],None,[],cfg.anchor_m)
    assert row['reconciliation_max_error']<1e-12
    assert row['totals']['branch_tail']>0.1
    assert row['totals']['blocked']==0
    assert row['totals']['unhit']==0
