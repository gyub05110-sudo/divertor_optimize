# Adaptive finalist peak comparison

The user authorized pushing the completed branch and proceeding with the next work described in the previous report: stop once rankings are proved, allocate work to unresolved comparisons, and avoid repeated independent replay of shared baselines. No additional approval is needed for this implementation.

## Decision

Keep production search and existing certificate semantics unchanged. Add an incremental certificate builder using the same directed-rounding operations and evidence schema. A shared comparison session owns one builder per kernel, repeatedly refines a participant in an unresolved comparison, and removes resolved comparisons from scheduling. Each returned final certificate is independently replayed once; all pairwise decisions are then recomputed from the verified immutable evidence. No caller-supplied success/cache flags are trusted.

Alternatives: repeating the existing certify call wastes earlier subdivisions; uniformly increasing budgets wastes effort on already separated candidates. Incremental shared scheduling is preferred because it retains all proved cells and gives actual early stopping.

The scheduling score is a heuristic: current interval width times the number of unresolved comparisons using that kernel. It affects runtime only, never correctness. Initial probes and worst-upper-cell subdivisions use the existing arithmetic. Ties are deterministic. Zero budget and overlapping intervals return unresolved. All source/resolution/problem identities must match within each pair; no cross-resolution intersection.

Generation has a shared wall-time and operation cap; deadline checks occur between atomic operations. Initialization, output creation, independent replay, and budget overrun are separately timed, with end-to-end time reported. Final replay is mandatory even if the generation budget expires. Runtime improvements are not claimed from fewer resolved decisions or looser evidence.

## Experiment and acceptance

Use exactly the six requalified 320 and 1280 kernels and five pair comparisons from the previous study. Run sequentially with the same 136-second total generation budget (two 8-second baseline budgets plus four 30-second candidate budgets previously). Preserve all prior evidence and input hashes. Compare resolved decisions, total time, generation operations, replay count, and interval brackets against the prior fixed-budget run. These remain finite frozen-kernel rankings, not continuous transport or global optimality proofs. Do not start a new optimization campaign or promote candidates automatically.
