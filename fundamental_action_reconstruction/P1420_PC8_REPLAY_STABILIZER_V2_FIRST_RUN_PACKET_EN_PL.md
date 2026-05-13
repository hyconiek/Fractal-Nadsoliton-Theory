# P1420 PC8 Replay Stabilizer v2 First Run Packet (EN/PL)

Status: `P1420_EXECUTED_PC8_RS2_RUN1_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

Run the first stabilized PC8 variant under unchanged strict thresholds.

## Execution

- Script: `p1420_pc8_replay_stabilizer_v2_first_run_checkpoint.py`
- Artifact: `generated/p1420_pc8_replay_stabilizer_v2_first_run_summary.json`

## Result

All preregistered thresholds passed in this run:

- boundary drift,
- selector margin,
- dual replay gap (at threshold boundary).

Verdict: `PASS_PC8_RS2_RUN1`.

## Lay explanation (PL)

Po ludzku: nowa wersja działa lepiej — testy stabilności i powtarzalności przeszły,
ale to nadal dopiero krok pośredni, nie końcowy dowód całej teorii.

## Recommendation

Execute strict selector-uniqueness discharge pre-check on this stabilized variant before any SM/GR promotion claims.
