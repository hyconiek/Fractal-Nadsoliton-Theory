# FIN unlimited local research run — P486, P485, P487, P475

This package runs exact local computations without an application-level time,
CPU, or memory limit. Operating-system and hardware limits still apply.

Execution is sequential:

1. `P486` — exact rational interval orientation/localization certificate.
2. `P485` — exact localized ideal-membership test for all five causal-axis
   consistency minors.
3. `P487` — structure-aware elimination after exact removal of `A,B,u` and
   selection of the positive standard branch.
4. `P475-unlimited` — original full 14-variable lexicographic elimination.

Sequential execution protects the Intel i3-10110U / 16 GB host from several
simultaneous exact-CAS workloads. There are no network calls, external models,
remote repositories, laboratory data, or external audits.

## Start

```bash
bash launch_fin_unlimited_research.sh
```

## Inspect without stopping

```bash
python3 fin_unlimited_pipeline_status.py
```

The master state is `FIN_Unlimited_Research_Pipeline_State.json`. Completion
creates `FIN_Unlimited_Research_Pipeline_DONE.json`. Each stage has a separate
log and status JSON. A stage failure is recorded; the pipeline then starts the
next independent stage.

Inside Codex, separate terminal calls may inhabit separate PID namespaces.
Consequently, a negative PID probe from the status helper is not by itself a
termination signal. The durable pipeline state, per-stage status, advancing
logs, and the unified-exec session recorded in
`FIN_Unlimited_Research_Exec_Session.txt` are authoritative.

Do not call an algebraic relation a minimal polynomial until its P473 factor
has been isolated, proved irreducible, and verified exactly. Do not call P485
proved unless all five exact remainders vanish and P486 has paid the nonzero
reference and orientation premises.
