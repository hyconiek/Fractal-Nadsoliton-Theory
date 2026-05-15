# P1800 S750 Strict BW Componentwise Divergence Trace Evidence Pack Packet (PL)

Status: `P1800_EXECUTED_STRICT_BW_COMPONENTWISE_DIVERGENCE_TRACE_EVIDENCE_PACK_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Uderzyć bezpośrednio w aktualny bottleneck BW:

1. ujednolicić format dowodu dla `EL_g - E_{μν}` i divergence trace,
2. wymusić komponentowe raportowanie `B1/B2/B3/C1/C2` na wspólnym freeze,
3. zablokować PASS bez jawnego zera residualu i jawnego zera divergence trace.

## Zakres strict-only

Brak legacy bridge i brak proxy. To jest wyłącznie strict-side evidence pack
na wejście do decyzji `G_BW`.

## Wymagany evidence pack

Każde wykonanie BW musi dostarczyć:

1. `freeze_id_common`, `background_family_id`, `index_convention_id`,
2. `component_residuals` dla `B1,B2,B3,C1,C2` z:
   - wartością symboliczną,
   - uproszczoną wartością końcową,
   - check digest,
3. `divergence_trace`:
   - `nabla_mu(E^{mu nu} - T^{mu nu})` w formie jawnej,
   - wynik uproszczenia końcowego,
   - check digest,
4. `classification` (`LOCAL/GLOBAL`, `NONPROXY`, `FULL_EXPORT`),
5. `verdict_candidate`: tylko `PASS_ZERO` lub `OPEN_OBSTRUCTION_WITH_TRACE`.

## Reguła werdyktu BW

`PASS_ZERO` jest dozwolony tylko gdy jednocześnie:

- wszystkie `B1,B2,B3,C1,C2` mają końcowe zero,
- divergence trace ma końcowe zero,
- brak freeze/index mismatch,
- pełny evidence pack ma check digests.

W każdym innym przypadku -> `OPEN_OBSTRUCTION_WITH_TRACE`.

## Co jest dowiedzione

1. BW ma teraz minimalny, kompletny format dowodowy na poziomie komponentów.
2. PASS bez full residual+divergence evidence jest formalnie zablokowany.
3. Pakiet jest kompatybilny z `P1797` (gate decision) i `P1799` (immutability audit).

## Co pozostaje OPEN

1. Realne podstawienie pełnych tensorów i uzyskanie zera dla wszystkich składowych.
2. Promocja do BRST/CUT po realnym `G_BW: PASS_ZERO`.

## Produkt

- Packet PL.
- Checkpoint JSON policy + evidence-pack template JSON.
