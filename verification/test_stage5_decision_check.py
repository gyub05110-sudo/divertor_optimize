import numpy as np
import pytest


def test_ordered_classification_reselects_wall_and_applies_precision_gate():
    from verification.stage5_decision_check import ordered_counts,compare
    n,rep=ordered_counts()
    assert n.shape==(462,6) and sum(rep)==46656
    a=np.arange(len(rep),dtype=float)
    # Identical winners against an alternative that is never lower.
    same=compare((5.,a+1),[(5.,a+1),(7.,a+3)],rep)
    assert same["estimate_pct"]==0 and same["ci95_pct"]==[0,0] and same["verdict_1pct"]=="supported"
    clear=compare((5.2,np.full(len(rep),5.2)),[(5.,np.full(len(rep),5.))],rep)
    assert clear["estimate_pct"]==pytest.approx(4.) and clear["verdict_1pct"]=="contradicted"
    wide=compare((5.,np.where(a%2,4.,6.)),[(5.,np.full(len(rep),5.))],rep)
    assert wide["verdict_1pct"]=="unresolved" and wide["uncertainty_pp"]==pytest.approx(20.)
