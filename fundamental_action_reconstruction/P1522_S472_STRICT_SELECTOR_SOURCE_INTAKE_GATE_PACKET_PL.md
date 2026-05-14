# P1522 S472 Strict Selector-Source Intake Gate Packet (No Legacy Bridge)

Status: `P1522_EXECUTED_STRICT_SELECTOR_SOURCE_INTAKE_GATE_NO_FALSE_CLOSURE`
As of: `2026-05-13`

## Cel

Wykonać kolejny uczciwy krok po `P1521/S471` na trasie:

```text
F_nadsoliton -> L_SM + L_GR
```

przez dodanie minimalnej bramki walidacyjnej `G_selector_intake^(strict)`.

Bramka nie domyka teorii; jedynie klasyfikuje, czy kandydat źródła selektora
może wejść do strict-core pipeline jako obiekt testowy.

## Kontrakt wejścia

Każdy kandydat źródła selektora musi dostarczyć:

1. `candidate_id`,
2. `provider_class`,
3. `symmetry_breaking_premise`,
4. `strict_provenance_trace`,
5. `noncyclic_anchor`.

Brak któregokolwiek pola => natychmiastowy `rejected`.

## Kontrakt wyjścia

`G_selector_intake^(strict)` zwraca wyłącznie:

1. `intake_status` in `{accepted_as_strict_source, rejected}`,
2. `reason_code`,
3. `provenance_hash`.

## Rygor i guardraile

1. Bez bridge do legacy, bez transferu legacy claims.
2. `QW-2191` pozostaje aktywny, dopóki nie ma zaakceptowanego strict-internal
   źródła selektora.
3. `QW-2381/2382/2383`: wymagany noncyclic anchor (`noncyclic_anchor=true`).
4. Bramka intake nie jest dowodem global closure; to tylko filtr wejścia.

## Kryterium PASS

PASS, gdy:

1. bramka odrzuca kandydatów niekompletnych lub cyklicznych,
2. bramka emituje audytowalny `reason_code`,
3. `provenance_hash` jest deterministyczny,
4. artefakt jawnie utrzymuje `qw2191_closed=false` bez zaakceptowanego źródła.

## Kryterium FAIL

FAIL, gdy:

1. kandydat bez strict proweniencji przechodzi jako accepted,
2. kandydat bez noncyclic anchor przechodzi jako accepted,
3. pojawia się claim selector closure bez nowego accepted strict source.
