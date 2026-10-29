# Independent IMP / OMP Optimization Implementation Plan

User approved the design regions and requested implementation with the historical optimizer on 2026-10-06. Work inline in the already isolated no-transport folder; preserve archive_math and historical runs.

Goal: independently normalize actual IMP and OMP sources and optimize the approved partial wall surfaces with the historical direct-mapping heat kernel, Sobol seeds, DE seeds and constraints.

Architecture: split the baseline polyline at selected region endpoints, retaining wall orientation and fixed fragments. Pass resulting whole split-segment targets to the preserved MultiSurfaceEvaluator. A separate evaluator replaces its initial legacy source profiles with actual midplane profiles and recomputes its baseline; no archived kernel modification. Each source has a separate search log, baseline, best candidate and refinement verdict.

- [ ] Add failing tests: source normalization/sign and OMP equivalence; split-wall selected lengths, fixed segments and source-independent design region membership; cached first hits equal complete candidate wall intersections.
- [ ] Implement `d3d_optimization/independent_sources.py`: reflected sampling adapter for IMP, physical source support, partial-segment splitting and per-source evaluator construction. Retain 4 B-spline coefficients per target, fixed endpoints, 120% length limit and shelf-only 1.5-degree floor.
- [ ] Implement `optimize_independent_sources.py`: IMP/OMP/both CLI, preflight, bounded smoke integration, production Sobol batches 1024/512/512 at 20/50/100 mm, two DE runs (seeds42/93, 64 members, 40 generations, no polish). Objective is worst whole-wall peak ratio across widths1/2/3/5/10mm. Count unique evaluated, geometrically valid and feasible cases separately. Never require or invent exactly10000 evaluations.
- [ ] Save actual timing, source normalization, approved selection policy, exact targets/fixed fragments, code provenance and coefficient vectors; reject overwrites. Smoke results are marked unqualified.
- [ ] Implement candidate re-evaluation at cubic320/2mm/161nodes, cubic640/1mm/321nodes and quintic640/1mm/321nodes, with independently reconstructed source and first-hit data. Do not promote failed candidates.
- [ ] Run focused new tests plus historical33 tests. Execute preflight and small search for each actual source, inspect failures, fixed-wall invariance, balance and direct whole-wall first-hit agreement. Run a finite microbenchmark separately; do not extrapolate old timing as a guarantee.
- [ ] Document launch and verification commands and historical timing evidence. Default full search remains runnable, but this implementation turn uses explicitly labelled small verification searches, not a claimed production optimum.

Historical evidence: screenshot run7233 unique candidates,4454 geometry-valid,2011 feasible,600.5646s. Expanded shelf-only run14702 unique candidates,13650 geometry-valid,2142 feasible,734.6987s (=12.24498min). New design has more than8 coefficients and two independent searches.
