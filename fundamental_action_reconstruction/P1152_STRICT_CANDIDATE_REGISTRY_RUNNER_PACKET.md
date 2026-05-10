# P1152 Strict Candidate Registry Runner Packet

Status: `P1152_EXECUTED_STRICT_CANDIDATE_REGISTRY_RUNNER_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wdrożyć następny uczciwy krok po `P1151`: porównywać wiele kandydatów przez ten
sam pipeline rygoru, zamiast oceniać je pojedynczo i niestandardowo.

## Professor-level decision

Dodaję runner rejestrowy, który:

1. bierze listę kandydatów JSON,
2. uruchamia dla każdego pełny `P1151` pipeline,
3. eksportuje tabelę `pass/fail`.

To daje porównywalność i ogranicza bias selekcyjny.

## Artifacts

- runner:
  `p1152_strict_candidate_registry_runner.py`
- sample registry:
  `generated/p1152_candidate_registry_example.json`
- sample failing candidate:
  `generated/p1152_candidate_input_fail_example.json`
- aggregate summary:
  `generated/p1152_strict_candidate_registry_runner_summary.json`

## Current run result

Dla przykładowego rejestru (1 pass + 1 fail):

```text
total=2, passed=1, failed=1
```

Runner zwraca kod błędu, gdy istnieje przynajmniej jeden fail — to celowe
zachowanie bezpieczeństwa rygoru.

## Honest boundary

`P1152` nie daje closure i nie rozwiązuje `QW-2191`; to narzędzie porównawcze
w rygorze metodologicznym.

## P1177 extension (`--rank-by-safe-margin`)

Dodano tryb `--rank-by-safe-margin`: po zakończeniu runu rejestru runner
automatycznie uruchamia `P1153` i dołącza shortlistę kandydatów o
najwyższym `safe_region_margin` (domyślnie top-3) do dalszych testów
fizycznych.

To jest kolejny uczciwy krok: najpierw pełny rygor pass/fail (`P1151`), potem
transparentna selekcja tylko wśród kandydatów, które przeszły pipeline i mają
obliczalny margines bezpiecznego regionu.
Dodatkowo dodano `--shortlist-output <path>`, aby zapisywać shortlistę do
osobnego artefaktu JSON dla downstream testów fizycznych bez ponownego
parsowania pełnego summary.


Tryb `--enforce-shortlist-consistency` uruchamia `P1177` po zbudowaniu shortlisty i
przerywa run przy niespójności, więc downstream testy fizyczne dostają wyłącznie
zweryfikowaną shortlistę.
