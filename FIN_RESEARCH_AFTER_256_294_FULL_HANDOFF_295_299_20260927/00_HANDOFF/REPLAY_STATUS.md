# FINAL REPLAY STATUS

Date: 2026-09-27

All five copied research scripts were executed again from `02_MACHINE_ARTIFACTS/` while assembling this handoff.

- 295 `complete_system_criterion_295.py` — exit 0; finite reversible-tape output reproduced.
- 296 `microscopic_periodic_bridge_296.py` — exit 0; periodic bridge output reproduced.
- 297 `fingerprint_optimal_probe_297.py` — `PASS`; `tau_star=0.5427059872630275`, worst-train Chernoff `0.006633619093729224`, held-out N8 `0.00871303017626598`.
- 298 `multitime_semigroup_fingerprint_298.py` — exit 0; N=6 multi-time/rate-drift output reproduced.
- 299 `memory_window_scaling_299.py` — `PASS`; `epsilon_mem(N=3)=0.049089770496779`, `epsilon_mem(N=8)=0.001893037885389`, drop factor `25.9317422412`.

The exact stdout snapshots from these final runs are stored as `02_MACHINE_ARTIFACTS/replay_stdout_295.txt` through `replay_stdout_299.txt`.
