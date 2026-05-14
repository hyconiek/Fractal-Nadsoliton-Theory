# P1580 / S530 Strict Kernel -> Full Chain -> Closure Gap Packet (PL)

Status: `P1580_EXECUTED_SEMIGLOBAL_T1579B_AND_W1578B_VALIDATION`
As of: `2026-05-14`

## Cel

Wykonać pełny tor strict-only bez bridge do legacy:

`K_strict -> współczynniki -> L_SM + L_GR -> równania ruchu`,

i dołączyć semiglobalną walidację:

1. `T1579B`: jednoznaczność `rho_gr` na pełnej domenie strict,
2. utrzymanie repliki `W1578B` poza punktami krytycznymi.

## Rygor

- Brak importu bridge legacy -> strict.
- Brak przenoszenia legacy physical-role claims.
- Priorytet strict-route: `F_Nadsoliton => L_SM + L_GR`.

## Artefakty wykonania

- `p1580_s530_strict_kernel_to_full_chain_and_closure_gap_checkpoint.py`
- `generated/p1580_s530_strict_chain_samples.csv`
- `generated/p1580_s530_strict_kernel_to_full_chain_and_closure_gap_summary.json`

## Wynik

Checkpoint eksportuje:

1. próbki `K_strict(d)` i pochodnej na gęstej domenie strict,
2. współczynniki strict `(c2, c4, cY)`,
3. szkic pełnego lagranżianu `L_SM + L_GR`,
4. szkic równań ruchu,
5. semiglobalny wynik `T1579B` i trwałość `W1578B` poza krytycznymi,
6. jawny status `strict_core_closure = OPEN` oraz listę brakujących obiektów:
   - strict internal selector source export,
   - symmetry-breaking witness,
   - strict selector uniqueness theorem,
   - full SM+GR global stability theorem.

## Następny uczciwy krok

Uruchomić `P1581`: zbudować `strict_internal_selector_source_export` oraz minimalny
`selector_symmetry_breaking_witness`, żeby przejść od kandydującego strict-chain
proof do formalnego postępu względem `QW-2191` bez sztucznego domknięcia.
