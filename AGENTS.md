# Project objective and working constraints

The user's objective is FAST divertor geometry optimization from a supplied baseline geometry and user-selected operating inputs. Hypothetical inputs are valid study inputs. Agreement with experimental measurements is NOT a project goal or a gate for implementation, optimization, or delivery.

- Keep the supplied real D3D baseline geometry fixed as the comparison reference. Do not change it to accommodate limitations of a deposition surrogate.
- Keep the requested operating inputs fixed during a run. Do not tune them merely to make the baseline pass.
- The primary optimization objective is to spread the arriving heat-load footprint and reduce its peak at fixed assigned power. Preserve toroidal winding in parallel connection-length integrals; do not substitute unexplained power attenuation or merely reward longer field lines for actual peak reduction.
- A missing-power failure of an approximate deposition kernel is a limitation of that kernel; do not describe it as evidence that the real D3D baseline is physically invalid.
- Diagnose omitted power routes, improve the transport/deposition representation, verify power accounting, arrival locations, heat-flux distribution and numerical convergence, then optimize with the same model.
- Build a reusable input-driven optimization tool, not a one-off model that works only at a convenient reference point. Define valid numerical input domains and report genuine numerical/geometry failures explicitly. Do not promise a feasible global optimum for every mathematically possible input.
- Prioritize runtime as an acceptance criterion. Record baseline evaluation and optimization runtimes; reuse fixed-field/path caches and deterministic evaluations.
- Establish search accuracy against the same model at higher numerical resolution and record cost at that accuracy. Use FIREFLY as a separate physics cross-check on fixed designs. Keep its outcomes out of the objective and report conflicting evidence with its operating conditions.
- Store new optimization studies, large data exports, code/input snapshots and accelerator caches on D:/DivertorOptimization rather than C:, unless the user requests another location. Use GPU acceleration when numerical parity and a useful measured speedup are established.
- Keep original models and saved studies reproducible. Prefer additive modules and new output directories; verify recorded hashes.
- Label reduced-model assumptions clearly. Conservation alone does not prove arrival locations are correct. Do not silently renormalize lost power, remove geometry constraints, cap widths, or equate numerical consistency with experimental validation.
- Do not add experimental calibration, measurement comparison, or full plasma-model validation to the required work unless the user requests them.

Superpowers skills: use for substantive implementation, planning, debugging, and significant verification. Skip for simple questions, file exploration, and explanations. Do not invoke using-superpowers automatically. Existing user authorization to execute the agreed workflow does not require repeated approval.
