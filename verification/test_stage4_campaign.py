from collections import Counter
from pathlib import Path

import pytest

from verification import stage3_campaign as common


def cases():
    return [dict(id=f"shot{i}", rank=i, status="ready", path=str(Path(f"{i}.parquet")),
                 sha256=f"{i:064x}", time_ms=1400.) for i in range(1,25)]


def test_schedule_has_independent_families_and_balanced_primary_order():
    from verification import stage4_campaign as campaign
    rows = cases()
    rows[4].update(status="missing_time",reason="no exact frame")
    rows[5].update(status="acquisition_failed",reason="pinned object unavailable")
    del rows[5]["path"]
    schedule = campaign.build_tasks(rows)
    assert all(j["equilibrium"] is None for j in schedule if j["input_rank"] == 6)
    assert len(schedule) == len({x["relative_output"] for x in schedule}) == 240
    assert Counter(x["family"] for x in schedule) == {
        "hotspot":72, "full":72, "floor-only":24, "floor-entry":24, "floor-shelf":24,"floor-connector":24}
    for seed in campaign.SEEDS:
        first = []
        for row in rows:
            group = [x for x in schedule if x["seed"] == seed and x["case_id"] == row["id"]]
            assert len(group) == (6 if seed == 30713 else 2)
            assert all(x["method"]=="screened" and x["budget_seconds"]==600 and
                       x["registration_status"]==row["status"] for x in group)
            first.append(group[0]["arm"])
        assert len(set(Counter(first).values())) == 1
    with pytest.raises(ValueError):
        campaign.build_tasks(rows[:-1])
    rows[2]["id"] = "../bad"
    with pytest.raises(ValueError):
        campaign.build_tasks(rows)


def test_worker_dispatch_accepts_new_registered_plan_and_recovers_excluded_job(tmp_path):
    from verification import stage4_campaign as campaign
    rows = cases()
    rows[0].update(status="missing_time",reason="no exact frame")
    jobs = campaign.build_tasks(rows)
    plan = dict(tasks=jobs,registration_sha256="aa"*32)
    lookup = lambda study: plan
    first = common.worker(tmp_path,0,plan_loader=lookup)
    assert first["status"] == "registration_excluded" and first["elapsed_s"] is None
    assert common.worker(tmp_path,0,plan_loader=lookup) == first
    assert not list((tmp_path/"attempts").rglob("launch.json"))
