# P1431 — QW-2191 (Release 8.1) Closure Status Clarification

Status: `P1431_QW2191_NOT_CLOSED_IN_STRICT_CORE_CLARIFIED`
As of: `2026-05-13`

## Why this packet exists

User feedback indicated confusion that the pipeline "closes QW-2191" and then appears to "close it again".
This packet clarifies the strict scientific status unambiguously.

## Professorial clarification (strict rigor)

`QW-2191` is **NOT closed in strict core** on current repo state.

What was closed in prior packets/checkpoints:

1. local replay stability edges,
2. artifact-level consistency checks,
3. partial certificate/discharge scaffolding,
4. non-strict or conditional closure edges.

What is **not** closed:

1. strict internal selector-source uniqueness discharge,
2. strict-core theorem-level selector closure,
3. strict projection-level closure for `F-Nadsoliton => L_SM + L_GR`.

## Single-source truth table

- `QW-2191 strict-core status`: `OPEN_UNIQUENESS_OBSTRUCTION`
- `global strict closure`: `NOT ACHIEVED`
- `allowed statement`: "local/staged closure artifacts exist"
- `forbidden statement`: "QW-2191 strict-core is solved"

## Interpretation rule for all future packets

If any checkpoint reports PASS, PASS may refer only to local contract scope.
PASS must never be interpreted as global strict-core closure unless explicit field says:

`strict_core_qw2191_closed = true`

and such field must be backed by a new exported strict theorem/witness.

## Next-step consequence

The correct next move is not "re-close QW-2191" rhetorically.
The correct next move is to construct either:

1. a strict internal selector source witness, or
2. a clearly marked non-strict selector premise branch.

Until then, QW-2191 remains open in strict core.

---

## PL (skrót)

Masz rację: nie wolno mówić, że QW-2191 jest zamknięte globalnie.
To, co było "PASS", dotyczyło lokalnych etapów (replay/stabilizacja), a nie domknięcia strict-core.

Wniosek: `QW-2191` nadal jest otwarte w strict-core i trzeba to jawnie utrzymywać
w każdym kolejnym pakiecie.
