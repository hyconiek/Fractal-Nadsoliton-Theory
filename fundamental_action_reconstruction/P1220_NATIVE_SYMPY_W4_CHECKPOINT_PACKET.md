# P1220 Native SymPy W4 Checkpoint Packet

Status: `P1220_EXECUTED_NATIVE_SYMPY_CHECKPOINT_WITH_HASH_CANONICALIZATION`
As of: `2026-05-11`

## Goal

Następny uczciwy krok: uruchomić natywny `P1205` z lokalnym SymPy i dopiąć
spójność hasha z `P1212`, aby kontrolowany `P1210` mógł być wykonany bez mocków.

## Professor-level decision

Aktualizuję `p1205_w4_sympy_cas_runner.py`:

1. `as_of -> 2026-05-11`,
2. hash liczony na kanonizacji zgodnej z `P1212`
   (`sort_keys=True`, `separators=(",",":")`, `ensure_ascii=False`).

## Honest boundary

To nadal checkpoint runtime, nie theorem-level discharge W4.
