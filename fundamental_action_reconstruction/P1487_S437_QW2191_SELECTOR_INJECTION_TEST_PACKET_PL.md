# P1487 — S4.37 QW-2191 Selector Injection Test (PL)

Status: `P1487_EXECUTED_QW2191_SELECTOR_INJECTION_TEST_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wstrzyknąć premise `SP1_SB_v1` do jawnego testu selektora i sprawdzić,
czy blokada unikalności QW-2191 znika **lokalnie** bez naruszenia policy gate,
w torze strict-only `F(nadsoliton) => L_SM + L_GR`.

## Decyzja profesorska

Wykonujemy test kontrastowy:

- baseline gap: `G0 = |W_SM - W_GR|`,
- injected gap: `G1 = |W_SM - W_GR - Delta_SB|`,
- kryterium postępu: `G1 < G0` oraz `|Delta_SB| <= safety_margin`.

Jeśli kryterium przechodzi, dostajemy lokalny dowód postępu selektorowego,
ale nie strict-final closure.
