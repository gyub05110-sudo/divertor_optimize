"""Independently audit recorded DIII-D base-budget trials. Does no thermal computation."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCES = ("OMP", "IMP")
SEEDS = (5161, 7753, 9923)
REVIEWED_DRIVER = "7940f20406c755b14ca3c52cc159348d9599b2a18f3536eedd727db7e1b01c28"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def identity(candidate):
    return hashlib.sha256(json.dumps(candidate, sort_keys=True, allow_nan=False).encode()).hexdigest()


def check(condition, message):
    if not condition:
        raise ValueError(message)


def events(rows, setup, policy, batch):
    yield "high_setup", None, setup["high_s"]
    if policy == "direct":
        for row in rows:
            yield "generation", row, row["generation_s"]
            yield "fine", row, row["fine"]["seconds"]
        return
    check(policy == "screen", f"unknown policy {policy}")
    yield "low_setup", None, setup["low_s"]
    for group in itertools.batched(rows, batch):
        for row in group:
            yield "generation", row, row["generation_s"]
            yield "coarse", row, row["coarse"]["seconds"]
        eligible = [row for row in group if row["coarse"]["eligible"]]
        if eligible:
            leader = min(eligible, key=lambda r: (r["coarse"]["peak"], r["i"]))
            yield "fine", leader, leader["fine"]["seconds"]


def replay_scores(rows, setup, policy, batch, budget):
    """Build an event timeline, truncate it, then rank only completed fine scores."""
    elapsed = 0.
    baseline = setup["baseline"]
    costs = dict(high_setup=0., low_setup=0., generation=0., coarse=0., fine=0.)
    result = dict(policy=policy, budget_s=budget, selected=None, peak=None, A2=None,
                  fine_ids=[], seen_coarse=[], eligible_fine=0, elapsed_last_completed_s=0.,
                  first_5pct_s=None, best_eligible_new=None)
    eligible_fine = []
    for kind, row, seconds in events(rows, setup, policy, batch):
        check(math.isfinite(seconds) and seconds >= 0., f"invalid {kind} cost")
        elapsed += seconds
        if elapsed > budget:
            break
        costs[kind] += seconds
        if kind in ("high_setup", "coarse", "fine"):
            result["elapsed_last_completed_s"] = elapsed
        if kind == "high_setup":
            result.update(peak=baseline["peak"], A2=baseline["A2"])
        if kind == "coarse":
            result["seen_coarse"].append(row["i"])
        if kind == "fine":
            result["fine_ids"].append(row["i"])
            if row["fine"]["eligible"]:
                eligible_fine.append((elapsed, row))
    result["eligible_fine"] = len(eligible_fine)
    if eligible_fine:
        _, winner = min(eligible_fine, key=lambda item: item[1]["fine"]["peak"])
        score = winner["fine"]
        result["best_eligible_new"] = dict(i=winner["i"], peak=score["peak"], A2=score["A2"])
        if score["peak"] < baseline["peak"]:
            result.update(selected=winner["i"], peak=score["peak"], A2=score["A2"])
        result["first_5pct_s"] = next((t for t, r in eligible_fine
                                      if r["fine"]["peak"] <= .95*baseline["peak"]), None)
    return result, costs


def historical_catalog(studies):
    """Keep source/radius identities as well as coefficient vectors across radii."""
    grid = studies / "nf-method-20261105/grid/design.json"
    entries = {source: [] for source in SOURCES}
    files = {}

    def records(path):
        if path.is_relative_to(ROOT):
            name = "repo/" + path.relative_to(ROOT).as_posix()
        else:
            name = "studies/" + path.relative_to(studies).as_posix()
        files[name] = sha(path)
        return read(path)

    def add(category, candidate):
        entries[candidate["source"]].append((category, candidate))

    for source, vectors in records(grid)["families"].items():
        check(len(vectors) == 28, f"{source} historical grid size")
        for v in vectors:
            add("28_point_grid", dict(source=source, radius_mm=20, vector=v))
    broad = studies / "imp-omp-global-20261029"
    extended = studies / "nf-method-20261105/campaign-2MW"
    for row in records(broad / "snapshot/campaign_warm_candidates.json")["candidates"]:
        add("warm_and_protected", row)
    for row in records(extended / "snapshot/campaign_protected_candidates.json")["candidates"]:
        add("warm_and_protected", row["candidate"])
    stages = studies / "omp-stage-qualification-20261104"
    for key, row in records(stages / "candidate_catalog.json").items():
        check(key == identity(row["candidate"]), "archived stage identity mismatch")
        add("published_stages", row["candidate"])
    paths = set((ROOT / "results").glob("*/*/R*/manifest.json"))
    for root in (broad, extended):
        paths.update(root.glob("*/R*/manifest.json"))
        paths.update(root.glob("*/R*/qualified/*/manifest.json"))
    paths.update(stages.glob("R*/qualified/*/manifest.json"))
    paths.update((studies / "omp-distribution-final-20261104").glob("qualified/*/manifest.json"))
    for path in sorted(paths):
        manifest = records(path)
        for entry in (manifest, manifest.get("search_best") or {}):
            if entry.get("candidate"):
                add("published_and_search_winners", entry["candidate"])
    for path in sorted((ROOT / "results/omp-tangent-curvature-20261029").glob("R*/manifest.json")):
        m = records(path)
        add("warm_and_protected", dict(source="OMP", radius_mm=m["radius_mm"], vector=m["coefficients_m"][:16]))
    m = records(ROOT / "results/regional-20261028/qualified/IMP/manifest.json")
    add("warm_and_protected", dict(source="IMP", radius_mm=100, vector=m["candidates"]["selected"][16:] + [.5]))
    return entries, files


def repository_vectors(files):
    """Use existing repository manifests and published grid anchors before storing any vectors."""
    candidates = {}
    for relative, expected in files.items():
        if not relative.startswith("repo/"):
            continue
        path = ROOT / relative.removeprefix("repo/")
        check(sha(path) == expected, f"changed published source {relative}")
        m = read(path)
        options = [entry["candidate"] for entry in (m, m.get("search_best") or {}) if entry.get("candidate")]
        if path.parent.parent.name == "omp-tangent-curvature-20261029":
            options.append(dict(source="OMP", radius_mm=m["radius_mm"], vector=m["coefficients_m"][:16]))
        if path.parent.name == "IMP" and path.parent.parent.name == "qualified":
            options.append(dict(source="IMP", radius_mm=100, vector=m["candidates"]["selected"][16:]+[.5]))
        for candidate in options:
            candidates[identity(candidate)] = candidate
    anchors = read(ROOT / "verification/fixtures/nf_search_anchors.json")
    for source, rows in anchors.items():
        vertices = [r["candidate"]["vector"] for r in rows]
        for i in range(7):
            for j in range(7-i):
                vector = [(a*(6-i-j)+i*b+j*c)/6 for a,b,c in zip(*vertices, strict=True)]
                candidate = dict(source=source, radius_mm=20, vector=vector)
                candidates[identity(candidate)] = candidate
    return candidates


def extract_history(studies):
    """Extract only comparison vectors missing from the repository and small baseline records."""
    entries, files = historical_catalog(studies)
    grid = studies / "nf-method-20261105/grid"
    broad = studies / "imp-omp-global-20261029"
    design = read(grid/"design.json")
    check(design["anchors_sha256"] == sha(ROOT/"verification/fixtures/nf_search_anchors.json"), "grid anchors changed")
    available = repository_vectors(files)
    memberships, extra = {}, {}
    for pairs in entries.values():
        for group, c in pairs:
            key = identity(c)
            memberships.setdefault(group, set()).add(key)
            if key not in available:
                extra[key] = c
    for name, expected in design["inputs"].items():
        check(sha(broad/"snapshot"/name) == expected, f"historical baseline input mismatch: {name}")
    original_geom = broad/"snapshot/input/D3D/optimized/d3d-divertor-floor-shelf/DIIID_geom.json"
    check(sha(original_geom) == read(broad/"snapshot/baseline_review/baseline_review.json")["baseline_geometry_sha256"],
          "archived geometry digest does not match the reviewed wall")
    baselines = {}
    for source in SOURCES:
        path = grid/f"{source}_chi0.3/result.json"
        old = read(path)
        baselines[source] = dict(file=f"nf-method-20261105/grid/{source}_chi0.3/result.json",
                                 file_sha256=sha(path), reference=old["reference"], baseline=old["baseline"])
    return dict(schema="nf-base-budget-history.v1", catalog_file_sha256=files, input_hashes=design["inputs"],
                anchors_sha256=design["anchors_sha256"], archived_grid_baselines=baselines,
                catalog_memberships={k:sorted(v) for k,v in memberships.items()}, additional_vectors=extra)


def load_history(studies):
    package = studies/"catalog.json"
    data = read(package) if package.is_file() else extract_history(studies)
    check(data["schema"] == "nf-base-budget-history.v1", "unknown history schema")
    check(sha(ROOT/"verification/fixtures/nf_search_anchors.json") == data["anchors_sha256"], "changed grid anchors")
    vectors = repository_vectors(data["catalog_file_sha256"])
    for key, c in data["additional_vectors"].items():
        check(identity(c) == key, "historical candidate identity mismatch")
        vectors[key] = c
    catalog = {src:[] for src in SOURCES}
    for group, keys in data["catalog_memberships"].items():
        for key in keys:
            c = vectors[key]
            catalog[c["source"]].append((group, c))
    return data, catalog


def freeze_evidence(study, folders):
    previous = study/"audit/audit.json"
    if previous.is_file():
        # Git checkouts and copied studies have new filesystem timestamps. Replay the
        # original observed order, binding each timestamp to the original bytes.
        stamp = read(previous)["freeze_order"]
    else:
        first = min((p for folder in folders for p in folder.glob("pair_*.json")),
                    key=lambda p:p.stat().st_mtime_ns)
        def capture(path):
            return dict(path=path.relative_to(study).as_posix(), sha256=sha(path), mtime_ns=path.stat().st_mtime_ns)
        stamp = dict(panels=[capture(p/"panel.json") for p in folders], first_saved_pair=capture(first))
    panels = stamp["panels"]
    first = stamp["first_saved_pair"]
    check({r["path"] for r in panels} == {(p/"panel.json").relative_to(study).as_posix() for p in folders},
          "incomplete panel timing evidence")
    check(first["path"] in {p.relative_to(study).as_posix() for d in folders for p in d.glob("pair_*.json")},
          "missing first pair")
    for item in [*panels, first]:
        check(sha(study/item["path"]) == item["sha256"], "freeze-order evidence hash mismatch")
    check(max(p["mtime_ns"] for p in panels) < first["mtime_ns"], "panel freeze must precede first saved pair")
    return stamp


def displacement(a, b, source):
    count = 8 if source == "IMP" else 16
    return max(abs(x-y) for x, y in zip(a[:count], b[:count], strict=True))


def novelty(rows, source, catalog):
    vectors = [r["candidate"]["vector"] for r in rows]
    separated = min(displacement(a, b, source) for a, b in itertools.combinations(vectors, 2))
    check(separated >= .001, f"{source}: proposals inside 1 mm duplicate tolerance")
    categories = {}
    for category in sorted({group for group, _ in catalog}):
        old = {identity(c): c["vector"] for group, c in catalog if group == category}
        exact_ids = sum(row["id"] in old for row in rows)
        exact_vectors = sum(v in old.values() for v in vectors)
        nearest = min((displacement(row["candidate"]["vector"], v, source), row["i"], key)
                      for row in rows for key, v in old.items())
        categories[category] = dict(compared_ids=len(old), exact_id_hits=exact_ids,
                                    exact_vector_hits=exact_vectors, min_coefficient_Linf_m=nearest[0],
                                    nearest_proposal=nearest[1], nearest_historical_id=nearest[2])
    return dict(min_within_panel_coefficient_Linf_m=separated, comparisons=categories)


def validate_measurement(state):
    check(isinstance(state["eligible"], bool), "invalid eligibility")
    if state["eligible"]:
        ledger = state["ledger"]
        check(not state["errors"] and all(math.isfinite(state[k]) and state[k] > 0 for k in ("peak", "A2")),
              "eligible score has errors or invalid peak/area")
        check(abs(ledger["prescribed"]-1.) <= 2e-12 and ledger["unhit"] <= .001
              and ledger["balance_error"] <= 2e-12
              and abs(ledger["prescribed"]-ledger["deposited"]-ledger["unhit"]) <= 2e-12,
              "eligible score has invalid power ledger")


def audit_panel(folder, setup, source, catalog):
    panel = read(folder / "panel.json")
    result = read(folder / "result.json")
    frozen = panel["rows"]
    check(len(frozen) == 12 and [r["i"] for r in frozen] == list(range(12)), f"{folder} panel size/order")
    check(len(list(folder.glob("pair_*.json"))) == 12, f"{folder} pair count")
    check(result["panel_sha256"] == sha(folder / "panel.json"), f"{folder} panel digest")
    check(panel["attempts"] == 12 + sum(panel["rejected"].values()) == frozen[-1]["cumulative_draws"]
          and all(a["cumulative_draws"] < b["cumulative_draws"] for a, b in itertools.pairwise(frozen)),
          f"{folder} draws/rejections")
    check(math.isclose(sum(r["generation_s"] for r in frozen), panel["generation_s"], abs_tol=1e-9)
          and panel["trailing_rejected_s"] >= 0., f"{folder} generation accounting")
    rows = []
    for item in frozen:
        row = read(folder / f"pair_{item['i']:02}.json")
        check(item == {k:row[k] for k in item} and identity(item["candidate"]) == item["id"],
              f"{folder} frozen candidate mismatch")
        c = item["candidate"]
        check(c["source"] == source and c["radius_mm"] == 20 and len(c["vector"]) == (9 if source == "IMP" else 16),
              f"{folder} source/radius/dimension")
        for level in ("coarse", "fine"):
            validate_measurement(row[level])
            check(math.isfinite(row[level]["seconds"]) and row[level]["seconds"] > 0., f"{folder} evaluation time")
        rows.append(row)
    measurement_hashes = {path.name: sha(path) for path in sorted(folder.glob("pair_*.json"))}
    fine = [r for r in rows if r["fine"]["eligible"]]
    winner = min(fine, key=lambda r:r["fine"]["peak"], default=None)
    baseline = setup["baseline"]["peak"]
    finite = min(baseline, winner["fine"]["peak"]) if winner else baseline
    check(result["finite_best_including_baseline"] == finite and result["candidates"] == len(rows)
          and result["fine_eligible"] == len(fine) and result["best_new"] == (winner["i"] if winner else None),
          f"{folder} exhaustive reference")
    summary = []
    completed_costs = {}
    for p in ("screen", "direct"):
        policy, costs = replay_scores(rows, setup, p, batch=4, budget=200.)
        completed_costs[p] = costs
        expected = result["policies"][p]
        for key, value in policy.items():
            equal = (math.isclose(value, expected[key], rel_tol=0., abs_tol=1e-8)
                     if isinstance(value, float) and expected[key] is not None else expected[key] == value)
            check(equal, f"{folder} {p}.{key}: {value} vs recorded {expected[key]}")
        summary.append(dict(source=source, seed=int(folder.name), policy=p, budget_s=200.,
                            charged_generation_completed_s=costs["generation"], coarse_n=len(policy["seen_coarse"]),
                            fine_n=len(policy["fine_ids"]), eligible_fine=policy["eligible_fine"],
                            selected=policy["selected"], peak=policy["peak"], A2=policy["A2"],
                            first_5pct_s=policy["first_5pct_s"], last_completed_s=policy["elapsed_last_completed_s"],
                            panel_best=finite, regret_pct=100*(policy["peak"]/finite-1) if policy["peak"] else None))
    return dict(panel_sha256=result["panel_sha256"], result_sha256=sha(folder/"result.json"),
                pair_sha256=measurement_hashes, draws=panel["attempts"], rejection_counts=panel["rejected"],
                generation_s=panel["generation_s"], trailing_after_last_accept_s=panel["trailing_rejected_s"],
                fine_eligible=len(fine), best_new=result["best_new"], novelty=novelty(rows, source, catalog),
                completed_event_costs_s=completed_costs, policy_records=result["policies"]), summary


def audit(study, studies):
    protocol = read(study / "protocol.json")
    folders = [study/s/str(n) for s in SOURCES for n in SEEDS]
    check(all((p/"result.json").is_file() for p in folders), "trial incomplete: all six results required")
    check(protocol["sources"] == list(SOURCES) and protocol["seeds"] == list(SEEDS)
          and [protocol[k] for k in ("size","batch","radius_mm","budget_s")] == [12,4,20,200],
          "unrecognized protocol")
    driver = ROOT / "verification/nf_base_budget.py"
    check(protocol["driver_sha256"] == sha(driver) == REVIEWED_DRIVER, "driver differs from reviewed generator")
    history, catalog = load_history(studies)
    check(protocol["inputs"] == history["input_hashes"], "different grid inputs")
    for p, digest in protocol["inputs"].items():
        check(sha(ROOT/p) == digest, f"different source input {p}")
    geometry_dir = ROOT/"input/D3D/optimized/d3d-divertor-floor-shelf"
    original_geom = geometry_dir/"DIIID_geom.json"
    geometry_provenance = read(geometry_dir/"manifest.json")["geometry"]
    wall_review = read(ROOT/"baseline_review/baseline_review.json")
    limiter = read(original_geom)["limiter"]
    check(sha(original_geom) == wall_review["baseline_geometry_sha256"] == geometry_provenance["sha256"]
          and wall_review["wall"] == [*limiter, limiter[0]], "supplied original limiter coordinates differ")
    stamp = freeze_evidence(study, folders)
    baselines, panels, summaries = {}, {}, []
    for source in SOURCES:
        setup = read(study/source/"setup.json")
        validate_measurement(setup["baseline"])
        old_grid = history["archived_grid_baselines"][source]
        previous = read(ROOT/"results/imp-omp-global-20261029"/source/"R050/manifest.json")
        check(setup["baseline"]["eligible"] and old_grid["reference"] == previous["reference"] == protocol["reference"],
              "baseline eligibility or operating-input mismatch")
        peak = setup["baseline"]["peak"]
        broad_peak = previous["baseline"]["peak_MW_m2"]
        check(math.isclose(peak, broad_peak, rel_tol=1e-12)
              and math.isclose(peak, old_grid["baseline"]["peak"], rel_tol=1e-12), "refined baseline peak mismatch")
        baselines[source] = dict(current_peak=peak, grid_peak=old_grid["baseline"]["peak"], broad_peak=broad_peak,
                                 low_setup_s=setup["low_s"], high_setup_s=setup["high_s"],
                                 setup_sha256=sha(study/source/"setup.json"))
        pooled = []
        for seed in SEEDS:
            folder = study/source/str(seed)
            panel, rows = audit_panel(folder, setup, source, catalog[source])
            panels[f"{source}/{seed}"] = panel
            summaries.extend(rows)
            pooled.extend((seed, row) for row in read(folder/"panel.json")["rows"])
        check(len({row["id"] for _, row in pooled}) == len(pooled), f"{source} repeated identities across seeds")
        baselines[source]["closest_between_seeds_Linf_m"] = min(
            displacement(a["candidate"]["vector"], b["candidate"]["vector"], source)
            for (i, a), (j, b) in itertools.combinations(pooled, 2) if i != j)
        parity = read(study/source/"parity.json")
        check(parity["id"] in {row["id"] for _, row in pooled} and 0 <= parity["relative_peak"] <= 1e-10,
              f"{source} missing or failed GPU/CPU peak check")
        baselines[source]["recorded_coarse_cpu_gpu_parity"] = parity
    return dict(driver_sha256=REVIEWED_DRIVER, protocol_sha256=sha(study/"protocol.json"), panel_hashes_verified=True,
                freeze_order=stamp,
                generation_review="Hash-pinned propose and panel draw seeded targets and project only against baseline "
                "geometry. All panels precede first saved candidate score. This is a source and artifact audit, "
                "not a runtime filesystem trace.", original_input_hashes=protocol["inputs"], baselines=baselines,
                original_geometry=dict(**geometry_provenance, original_limiter_vertices=len(limiter),
                                       exact_coordinate_match_with_added_closure=True),
                distance_definition="Maximum absolute difference in metre-valued coefficients; the dimensionless "
                "IMP junction fraction is included in exact vector equality and excluded from metre distances.",
                catalog_scope="Saved 28-point grids, named qualified/search winners, published OMP stage vectors, "
                "campaign warm/protected vectors. Excludes all-time proposal logs.",
                catalog_file_sha256=history["catalog_file_sha256"],
                eligibility_scope="Recorded eligibility used in replay; eligible power ledger, errors and positive "
                "peak/area checked. Curves and thermal fields not independently recomputed.", panels=panels), summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, help="new trial directory containing protocol.json")
    parser.add_argument("--studies", type=Path, help="archive study root, or portable directory containing catalog.json")
    parser.add_argument("--export-history", type=Path, help="create a compact catalog from the archived --studies")
    parser.add_argument("--write", action="store_true", help="publish audit/audit.json and summary.csv once")
    args = parser.parse_args()
    if args.export_history:
        if not args.studies or args.write or args.study:
            parser.error("--export-history requires --studies alone")
        history = extract_history(args.studies.resolve())
        args.export_history.mkdir(parents=True, exist_ok=True)
        with (args.export_history/"catalog.json").open("x", encoding="utf8", newline="\n") as stream:
            json.dump(history, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        return
    if not args.study:
        parser.error("--study is required for the audit")
    study = args.study.resolve()
    report, rows = audit(study, args.studies or study.parent)
    out = study / "audit"
    if not args.write and (out/"audit.json").is_file():
        check(report == read(out/"audit.json"), "replay differs from the published audit or paired-file hashes")
        with (out/"summary.csv").open(newline="", encoding="utf8") as stream:
            published = list(csv.DictReader(stream))
        expected = [{k:"" if v is None else str(v) for k,v in r.items()} for r in rows]
        check(published == expected, "recomputed policy summary differs from published CSV")
    if args.write:
        if out.exists():
            raise FileExistsError(out)
        with tempfile.TemporaryDirectory(dir=study) as tmp:
            ready = Path(tmp)/"audit"
            ready.mkdir()
            (ready/"audit.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf8")
            with (ready/"summary.csv").open("w", newline="", encoding="utf8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            ready.rename(out)
    print(f"Verified {len(report['panels'])} panels and {len(rows)} matched policy results.")


if __name__ == "__main__":
    main()
