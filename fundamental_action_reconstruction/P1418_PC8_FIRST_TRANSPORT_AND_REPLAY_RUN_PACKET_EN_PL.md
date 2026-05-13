# P1418 PC8 First Transport and Replay Run Packet (EN/PL)

Status: `P1418_EXECUTED_PC8_RUN1_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

Execute first PC8 run with unchanged strict thresholds from design freeze.

## Execution

- Script: `p1418_pc8_first_transport_and_replay_run_checkpoint.py`
- Summary: `generated/p1418_pc8_first_transport_and_replay_run_summary.json`
- Obstruction: `generated/p1418_pc8_replay_gap_obstruction_v1.json`

## Result

- Boundary drift passed,
- selector margin passed,
- dual replay gap failed (`FAIL_PC8_REPLAY_GAP`).

Interpretation:

PC8 improves selector separation but still lacks replay stability needed for strict promotion.

## Lay explanation (PL)

Dla laika: model już lepiej „odróżnia” najlepszą odpowiedź i jest stabilny na zmianę warunków,
ale dwie niezależne powtórki nadal różnią się odrobinę za bardzo.
Czyli jesteśmy bliżej celu, ale jeszcze nie mamy wystarczającej powtarzalności.

## Recommendation

Build `PC8_replay_stabilizer_v2_noncyclic` and rerun P1418-equivalent checks with the same thresholds.
