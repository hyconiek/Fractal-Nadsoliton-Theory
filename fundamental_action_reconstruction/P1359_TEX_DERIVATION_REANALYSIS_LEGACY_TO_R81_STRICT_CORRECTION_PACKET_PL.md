# P1359 TeX Derivation Reanalysis: Legacy -> R8.1 Strict Correction Packet (PL)

Status: `P1359_EXECUTED_TEX_REANALYSIS_AND_CORRECTION_NO_FALSE_PASS`
As of: `2026-05-12`
Inputs:
- `TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex`
- `TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.tex`
- `RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md`

## Cel

Uzupełnić brak: jawnie przeanalizować wyprowadzenia opisane w TeX i skorygować je do stanu R8.1 strict.

## Co pokazuje TeX legacy 4.4

W legacy TeX są bezpośrednie, mocne wyprowadzenia typu:

1. `sin^2(theta_W)=alpha_geo/12` jako centralna formuła,
2. jawne wartości `g1,g2,g3` jako „derived from kernel”,
3. liczne twierdzenia o emergencji `SU(3)xSU(2)xU(1)` i o zgodności grawitacyjnej,
4. narracje o „UV-completion of gravity”.

To jest historycznie ważne, ale metodologicznie jest to warstwa legacy, nie automatycznie strict-R8.1.

## Co pokazuje TeX strict R8

W strict TeX (R8) wielokrotnie zapisano hard limits:

1. brak automatycznego transferu legacy-claims na strict kernel,
2. brak claimu pełnego kernel-alone/global discharge `QW-2191` w wielu miejscach dokumentu,
3. rozróżnienie: scaffold/packaging vs pełny theorem-level dowód fizyczny,
4. otwarte granice dla pełnej ToE closure i pełnej fizycznej identyfikacji.

## Korekta wcześniejszych wyprowadzeń (z TeX legacy) po R8.1

| Klasa wywodu w TeX 4.4 | Korekta po R8.1 strict |
|---|---|
| `sin^2(theta_W)=alpha_geo/12` jako finalny dowód | traktować jako legacy-identity; w strict wymaga osobnego successor/provenance map |
| `g1,g2,g3` „z kernela” | utrzymać jako candidate lane, ale nie jako finalną zgodność bez residual tables i provenance lock |
| `SU(3)xSU(2)xU(1)` „emerges fully” | w strict: mocny scaffold, ale nie pełna fizyczna unikalność końcowa |
| „gravity UV completion established” | w strict R8.1: nie wolno promować do finalnego claimu bez współczesnego strict residual/audit bridge |

## Odpowiedź na Twoje pytanie (po analizie TeX)

Czy TeX-owe wyprowadzenia już dowodzą, że z kernela wychodzą znane wartości fizyczne?

**Nie w aktualnym strict sensie dowodowym.**

Po R8.1 one stanowią:

- mocny materiał startowy + hipotezy liczbowe,
- ale wymagają ponownego strict pipeline replay (`P1358` pokazał pierwszy realny FAIL na niekalibrowanej mapie),
- więc nie można ich bezpośrednio traktować jako domknięty dowód ToE.

## Decyzja profesorska

Najuczciwszy następny krok po tej analizie TeX:

1. zbudować `P1360` = TeX-to-Strict Claim Compiler,
2. dla każdej legacy formuły automatycznie przypisać status: `strict_verified / strict_candidate / legacy_only`,
3. automatycznie generować residual-obligation listę dla pozycji `strict_candidate`,
4. publikować tylko te tezy, które mają `strict_verified` + residual pass.

## Dla laika

Stare dokumenty TeX zawierają dużo trafnych intuicji i dobrych pomysłów.
Ale po nowym rygorze R8.1 trzeba każdą starą tezę przepuścić przez nowy „egzamin jakości”.
Dopiero te, które go przejdą, mogą być nazywane prawdziwie potwierdzonym wynikiem fizycznym.
