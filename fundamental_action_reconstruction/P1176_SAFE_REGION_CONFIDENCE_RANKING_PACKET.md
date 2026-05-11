# P1176 Safe-Region Confidence Ranking Packet

Status: `P1176_EXECUTED_SAFE_REGION_CONFIDENCE_RANKING_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać następny krok po `P1175`: rozszerzyć ranking o miarę pewności
bezpiecznego regionu, a nie tylko pass/fail.

## Professor-level decision

Rozszerzam `P1153` o metrykę:

```text
safe_region_margin = min distance do granic safe_bounds
```

Interpretacja:

- większy `safe_region_margin` = większy zapas bezpieczeństwa,
- `-1.0` = kandydat poza regionem,
- `None` = brak danych hint.

Sortowanie rankingu uwzględnia teraz ten margines po `quality_score`.

## Result

Zaktualizowany ranking eksportuje pole `safe_region_margin` dla każdego
kandydata i lepiej rozróżnia kandydatów "ledwo w regionie" od stabilnych w
rdzeniu regionu.

## Artifacts

- updated ranking script:
  `p1153_strict_candidate_quality_ranking.py`
- updated summary:
  `generated/p1153_strict_candidate_quality_ranking_summary.json`

## Honest boundary

`P1176` to ulepszenie decyzji operacyjnych, nie claim closure ani discharge
`QW-2191`.
