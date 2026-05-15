# P1786 S736 Strict Current Priority Nonproxy Closure Gate Matrix Packet (PL)

## Cel

Ustalić jeden jawny, strict-only gate-matrix dla aktualnego bottlenecka:
`E_A^μ`, `E_H`, `EL_g`, boundary-control, `H1` 4D weak-form,
Bianchi/Ward, BRST, Cutkosky.

## Zasada strict-only

- Zero legacy bridge.
- Zero proxy shortcut.
- Zero PASS bez jawnego witness/residual.

## Gate matrix (stan na teraz)

1. `G1_EA_NONPROXY_EXPLICIT_EXPORT` -> `OPEN`
   - wymagane: indeksowo-komponentowe `E_A^μ` na wspólnej rodzinie teł.
2. `G2_EH_NONPROXY_EXPLICIT_EXPORT` -> `OPEN`
   - wymagane: indeksowo-komponentowe `E_H` na tej samej rodzinie teł.
3. `G3_ELG_NONPROXY_EXPLICIT_EXPORT` -> `OPEN`
   - wymagane: komponentowy `EL_g^{μν}` i jawny residual `EL_g-E_{μν}`.
4. `G4_BOUNDARY_TERM_CONTROL` -> `PARTIAL_OPEN`
   - wymagane: domknięcie boundary-control dla tej samej konfiguracji testowej.
5. `G5_H1_4D_WEAK_FORM` -> `OPEN`
   - wymagane: H1 witness `δE_A^μ/δH - δE_H/δA_μ` (4D, nonproxy).
6. `G6_BIANCHI_WARD` -> `OPEN`
   - wymagane: divergence trace + zgodność z tożsamościami.
7. `G7_BRST_NILPOTENCY` -> `OPEN`
   - wymagane: jawny nilpotency witness na aktualnym bundle.
8. `G8_CUTKOSKY_UNITARITY` -> `OPEN`
   - wymagane: jawny cut-compatibility witness.

## Co jest dowiedzione (lokalnie)

- Istnieją jawne definicje nonproxy `E_A^μ`, `E_H`, `EL_g^{μν}` jako targety wariacyjne.
- Istnieje strict-only sygnalizacja, że theorem-level QG closure jest nadal `OPEN`.

## Co jest OPEN (globalnie)

- Brakuje pełnych komponentowych eksportów i residual zero witness dla `E_A^μ/E_H/EL_g`.
- Brakuje H1, Bianchi/Ward, BRST i Cutkosky witness na poziomie theorem-gate.

## Ryzyko false-pass

- Oznaczenie `FULL EXPORT` przy samych template'ach.
- Oznaczenie `PASS` bez residual vector / witness ledger.
- Mieszanie stanów `LOCAL/REDUCED` z `GLOBAL/NONPROXY`.

## Następny uczciwy krok

Wykonać `G1+G2` razem na tej samej rodzinie teł (jeden run-contract),
a następnie uruchomić `G5` (H1) bez claimu `PASS` do czasu jawnego residual=0.

## Objaśnienie dla laika

To jest lista bramek jakości: dopóki każda bramka nie ma twardego dowodu,
teoria pozostaje uczciwie oznaczona jako "jeszcze otwarta".
