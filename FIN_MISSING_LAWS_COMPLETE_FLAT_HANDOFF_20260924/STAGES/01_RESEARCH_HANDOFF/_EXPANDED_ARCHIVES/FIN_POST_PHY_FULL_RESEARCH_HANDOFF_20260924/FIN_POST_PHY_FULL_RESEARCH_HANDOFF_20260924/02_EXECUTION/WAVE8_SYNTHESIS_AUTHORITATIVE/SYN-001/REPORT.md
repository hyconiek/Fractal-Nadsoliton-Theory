# SYN-001 — adversarial proof and scope consolidation

Status: **PARTIAL_REVIEW_READY__EARLY_WAVE_ARTIFACT_PROVENANCE_GAP**.

The current runtime contains full task packages for Waves 4–7, the Wave-3 summary, the original programme/specification, and the authoritative repository sources. It does **not** contain the full per-task Wave-0/1/2 artifact directories produced earlier in the conversation. Those tasks were previously counted as completed, but their byte-level manifests and replay scripts are not available here for an independent SYN-001 recheck.

Available execution material: FIN_POST_PHY_EXEC_20260924_WAVE3, FIN_POST_PHY_EXEC_20260924_WAVE4, FIN_POST_PHY_EXEC_20260924_WAVE5, FIN_POST_PHY_EXEC_20260924_WAVE6_X, FIN_POST_PHY_EXEC_20260924_WAVE7.

## Adversarial checks on currently available artifacts

- no new scientific GitHub commit intervenes: main remains `d4c4a0ac7933ff53a958976ce6ddbde5cd514017`, child of baseline `01c89aa3c24ce5840ef86e950530b6ae8caafd24`;
- new class-X results preserve source-vs-premise boundaries;
- REF-006 explicitly separates fixed-band local convergence from tail/fiber nonuniqueness and from fine modes;
- GEO-002 labels its stationary mixture non-ergodic instead of misusing it as an ergodic theorem;
- DYN-004 uses the declared dissipative kinetic law and does not add inertia/drive after the no-drift result;
- GATE-002 remains closed;
- external empirical tasks remain NOT_READY/NO_FEASIBLE_PILOT rather than simulation-as-data.

## Quarantine

Because Wave-0/1/2 task bytes are unavailable in this runtime, SYN-001 cannot honestly issue a blanket PASS over all 42 dependency paths. Claims depending only on currently auditable repository sources and Waves 3–7 may retain their scoped outcomes; a *global* review-ready promotion is quarantined until the early execution packages are restored or regenerated from their exact task specifications.

This is a successful fail-closed audit outcome, not evidence against the mathematics of those earlier tasks.
