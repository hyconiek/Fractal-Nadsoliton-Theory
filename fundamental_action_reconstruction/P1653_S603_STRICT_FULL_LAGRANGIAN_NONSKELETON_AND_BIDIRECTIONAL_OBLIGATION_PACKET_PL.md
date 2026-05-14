# P1653 / S603 — Strict full Lagrangian (non-skeleton) and bidirectional obligation packet

## Cel
Wykonać kolejny uczciwy krok bez bridge do legacy:
1. wyeksportować **pełną (nieszkieletową) postać** gęstości lagranżianu strict,
2. pokazać jawny tor `K_strict -> współczynniki -> L_total -> EOM`,
3. zapisać obowiązki dla toru odwrotnego `EOM -> L_total -> współczynniki -> K_strict` bez false-pass.

## Zakres strict-only
- `strict_only = true`,
- `legacy_bridge_used = false`,
- bez translacji roszczeń legacy (`alpha_geo/beta_tors`) na strict-kernel.

## Wynik merytoryczny
- Używamy jawnego Lagrangianu z sektorami:
  - skalar strict (`phi`) z potencjałem wielomianowym,
  - pełny sektor cechowania SM (`SU(3)xSU(2)xU(1)`),
  - fermiony SM z Yukawami,
  - Higgs,
  - sektor grawitacyjny z członami wyższej krzywizny,
  - sprzężenia mieszane (`phi-H`, `phi-R`, `H-R`).
- Definiujemy mapę obowiązków EL/Helmholtza dla toru odwrotnego jako jawny backlog theorem-level.

## Wyjście
- `generated/p1653_s603_strict_full_lagrangian_nonskeleton_and_bidirectional_obligation_summary.json`
