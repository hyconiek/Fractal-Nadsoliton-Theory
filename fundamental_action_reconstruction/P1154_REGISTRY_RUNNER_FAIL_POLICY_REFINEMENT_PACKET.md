# P1154 Registry Runner Fail Policy Refinement Packet

Status: `P1154_EXECUTED_REGISTRY_RUNNER_FAIL_POLICY_REFINEMENT_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Usprawnić rygor operacyjny: odróżnić tryb CI (fail-fast) od trybu badawczego
(agreguj wyniki i kontynuuj ranking), bez rozluźniania guardraili.

## Professor-level decision

`P1152` dostaje opcję:

```text
--allow-failures
```

Semantyka:

1. domyślnie pozostaje fail-fast (`exit 1` jeśli jakikolwiek kandydat fail),
2. z `--allow-failures` runner kończy się `exit 0`, ale nadal raportuje
   `failed > 0` w summary JSON,
3. dzięki temu można uczciwie wykonać `P1153` ranking także dla mieszanych
   rejestrów (pass+fail) bez ukrywania porażek.

## Outcome

Refinement zachowuje rygor (porażki są jawne), a jednocześnie usuwa
techniczny blocker sekwencji `registry -> ranking`.

## Honest boundary

To zmiana czysto metodologiczno-operacyjna. Nie jest twierdzeniem fizycznym,
nie rozwiązuje `QW-2191` i nie daje closure.
