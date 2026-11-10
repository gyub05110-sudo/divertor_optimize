import pytest
from verification import stage5_campaign as study


def make_roles():
    labels = (("baseline", None), *((a,30713) for a in study.ARMS))
    records = []
    for j, case in enumerate(study.CASES):
        roles = [dict(case_id=case,role=role,optimizer_seed=seed,wall_id=f"{i:064x}")
                 for i,(role,seed) in enumerate(labels)]
        if case==study.EXCEPTION:
            roles.extend(dict(case_id=case,role=a,optimizer_seed=30757,wall_id=f"{10+i:064x}")
                         for i,a in enumerate(("hotspot","full")))
        records.extend(roles)
    return records


def test_schedule_covers_exact_matrix_shares_only_within_case():
    rows = make_roles()
    rows[2]["wall_id"] = rows[1]["wall_id"]
    jobs = study.jobs(rows)
    assert len(jobs) == 50*6
    assert len({(c,w,s) for c,w,s in jobs}) == len(jobs)
    assert set(s for _,_,s in jobs)==set(study.SEEDS)
    assert {c for c,_,_ in jobs}==set(study.CASES)
    with pytest.raises(ValueError):
        study.jobs(rows[:-1])
    with pytest.raises(ValueError):
        study.jobs(rows + [rows[0]])


def test_panel_check_rejects_swapped_geometry():
    rows=make_roles()
    panel={"cases":[dict(case_id=c,roles=[dict(role=r["role"],seed=r["optimizer_seed"],wall_id=r["wall_id"])
                                                  for r in rows if r["case_id"]==c]) for c in study.CASES]}
    study.verify_panel(rows,panel)
    panel["cases"][0]["roles"][0]["wall_id"]="f"*64
    with pytest.raises(ValueError):
        study.verify_panel(rows,panel)
