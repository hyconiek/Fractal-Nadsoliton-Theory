# P1765 / S715 — Jawny nonproxy eksport metryczny `EL_g^{μν}`

Status: `P1765_S715_EXPLICIT_NONPROXY_METRIC_ELG_OPERATOR_EXPORT_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Wykonano następny bottleneck po `P1764`: dodano jawny, strict-only eksport
metrycznego operatora wariacyjnego:

```text
EL_g^{μν} := (2/sqrt(-g)) δS/δg_{μν}
```

na torze:

```text
K_strict -> coefficients -> full non-skeleton L_total -> metric nonproxy EL_g export
```

z pełnym zakotwiczeniem sektorowym (`gravity`, `gauge`, `fermion`, `Higgs`,
`scalar_phi`, `mix`).

## Co zostało dowiedzione

Na poziomie operatorowym (kowariantnym) jawnie wyodrębniono strukturę:

1. część grawitacyjna: `G^{μν}`, `Λ g^{μν}`,
2. poprawki krzywiznowe: `H_R2^{μν}`, `H_Ricci2^{μν}`, `H_Riemann2^{μν}`,
3. tensory energii-pędu sektorów materii i miksów,
4. wkłady kontrterminowe `T_CT^{μν}`.

To formalnie odblokowuje wejście do residual testu:
`EL_g - E_{μν}` na bazie `B1/B2/B3/C1/C2`.

## Co nadal jest OPEN

1. Brak komponentowego (indeksowego) policzenia residualu `EL_g - E_{μν}`.
2. Brak jawnego divergence trace dla bramki Bianchi/Ward.
3. Brak theorem-level witnessów BRST/Cutkosky.
4. Brak globalnego Helmholtz integrability closure.

## Ryzyka false-pass

1. Mylenie eksportu operatorowego z policzonym residualem komponentowym.
2. Nadanie statusu `PASS` bez jawnego wektora residuali i śladu obstrukcji.
3. Traktowanie gotowości weak-form jako strict-core closure.

W tym pakiecie nie wydano `PASS_ZERO`.

## Następny uczciwy krok

1. Policz `EL_g - E_{μν}` komponentowo na `B1/B2/B3/C1/C2`.
2. Opublikuj tylko: `PASS_ZERO` albo `OBSTRUCTION` + trace składników.
3. Na tych samych danych dołącz divergence check (`∇_μ E^{μν}`) dla
   Bianchi/Ward.
4. Dopiero potem aktualizuj bramki BRST/Cutkosky.

## Dla laika

To krok podobny do zapisania pełnego równania sił w fizyce: wiemy już, jakie
człony muszą wystąpić. Ale prawdziwy test dopiero przed nami: trzeba policzyć
każdy składnik i sprawdzić, czy równowaga naprawdę wychodzi dokładnie.
