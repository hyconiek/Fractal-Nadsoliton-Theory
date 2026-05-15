# P1769 / S719 — Pierwsza próba wykonania G_BW (componentwise)

Status: `P1769_S719_GBW_FIRST_COMPONENTWISE_EXECUTION_OBSTRUCTION_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Wykonano pierwszą pełną próbę uruchomienia bramki `G_BW` na przygotowanym
input-packu `P1768`, zgodnie z sekwencją `P1767` i bez zmiany klasy tła.

Werdykt tej próby:

`OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Co zostało dowiedzione

1. Pipeline strict-only przeszedł z etapu „gotowości” do realnego wykonania
   komponentowego testu Bianchi/Ward.
2. Obstrukcja została zlokalizowana jawnie w ciężkich członach tensorowych:
   - `∇_μ(H_R2^{μν})`,
   - `∇_μ(H_Ricci2^{μν})`,
   - `∇_μ(H_Riemann2^{μν})`,
   - oraz w nierozwiniętej bazie `T_CT^{μν}`.
3. Nie wydano żadnego `PASS_ZERO`.

## Co nadal jest OPEN

1. Domknięcie komponentowego residualu `EL_g-E_{μν}` dla `B1/B2/B3/C1/C2`.
2. Domknięcie divergence trace do identycznego zera.
3. Odblokowanie `G_BRST` i `G_CUT` (nadal formalnie zablokowane przez `G_BW`).

## Ryzyka false-pass

1. Zastąpienie nierozwiniętych tensorów deklaracją „prawdopodobnie zero”.
2. Lokalny PASS w reduced-subsectorze mylony z globalnym nonproxy closure.
3. Ukrycie wkładu kontrterminów (`T_CT`) bez jawnego śladu składnikowego.

## Następny uczciwy krok

1. Rozwinąć komponentowo brakujące człony `H_R2/H_Ricci2/H_Riemann2`.
2. Wyeksportować jawny tensorowy basis dla `T_CT^{μν}` w tej samej konwencji.
3. Powtórzyć `G_BW` na tym samym tle z dokładnie tym samym kontraktem wyjścia:
   `PASS_ZERO` albo `OBSTRUCTION_WITH_DIVERGENCE_TRACE`.

## Dla laika

To był prawdziwy „test silnika” teorii: uruchomiliśmy procedurę i wiemy już,
gdzie dokładnie jest blokada obliczeń. To postęp, bo problem jest teraz
konkretny i można go systematycznie usunąć bez zgadywania.
