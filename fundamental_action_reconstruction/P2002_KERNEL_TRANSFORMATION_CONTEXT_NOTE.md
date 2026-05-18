# P2002 Kernel Transformation Context Note (Legacy provenance, strict usage discipline)

Ten krótki note powstał po przeglądzie `DIAGRAMS_KERNEL_TRANSFORMATION.md`.

## Co wynika z diagramu historycznie

W diagramie legacy kernel powstaje jako efektywna kompresja trasy
`K_total -> K(d)=alpha*cos(omega d + phi)/(1+beta d)` z naciskiem na mechanizmy
rezonansowo-topologiczne i hiperboliczne tłumienie.

## Jak to używać dziś bez naruszenia guardrails

1. Traktować ten materiał jako **kontekst historyczny/proweniencyjny** dla legacy lane.
2. W strict lane (P1997+) NIE przenosić automatycznie ról legacy do `K_strict_gate`.
3. Utrzymać rozdział: legacy-provenance context vs strict operational witness pipeline.

## Znaczenie dla P2002

P2002 nie używa legacy równań do obliczeń; używa tylko strict-side artifactów
(P2000/P2001) do klasyfikacji trendu `Delta_opt`.
