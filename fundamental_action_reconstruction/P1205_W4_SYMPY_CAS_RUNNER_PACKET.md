# P1205 W4 SymPy CAS Runner Packet

Status: `P1205_PREPARED_EXTERNAL_SYMPY_RUNNER_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Dostarczyć skrypt do uruchomienia w środowisku z CAS (`sympy`), który generuje
plik wynikowy gotowy do wrzucenia na GitHub jako artefakt.

## Script

`p1205_w4_sympy_cas_runner.py`

Domyślny plik wynikowy:

`fundamental_action_reconstruction/generated/p1205_w4_sympy_cas_runner_summary.json`

## Run command

```bash
python3 fundamental_action_reconstruction/p1205_w4_sympy_cas_runner.py
```

Opcjonalnie własna ścieżka:

```bash
python3 fundamental_action_reconstruction/p1205_w4_sympy_cas_runner.py \
  --out /tmp/p1205_w4_sympy_cas_runner_summary.json
```

## Honest boundary

Skrypt generuje artefakt pre-discharge (expand/collect/simplify/hash) i jawnie
nie wykonuje discharge W4 ani nie zgłasza closure claim.
