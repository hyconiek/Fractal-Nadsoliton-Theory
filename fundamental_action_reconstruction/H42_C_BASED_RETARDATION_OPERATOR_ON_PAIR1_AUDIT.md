# H42 c-Based Retardation Operator On Pair1 Audit

Status: `PASS_PARTIAL_C_BASED_RETARDATION_SPLIT_ESTABLISHED`
Date: `2026-03-07`

## Goal

Zdefiniowac minimalny operator retardacji oparty o jawne `c` na
`pair1 = (c_1, s_1)` i sprawdzic dwa przypadki:

1. bez importowanego `psi0`,
2. z importowanym `psi0`.

Pytanie robocze brzmi:
- czy samo `c` daje selector mechanism,
- czy tylko kontroluje szybkosc/faze sprzezenia,
- i czy nowy efekt pojawia sie dopiero po imporcie anchoru `psi0`.

## Pair carrier

Pracujemy na:

- `V_1 = span{c_1, s_1}`

z baza:

- `(c_1, s_1)`.

## Minimal c-based retardation phase

Niech:

- `phi_ret(L;c) = omega * L / c`

bedzie minimalna faza retardacji dla drogi o dlugosci `L`.

## Case A: no psi0 imported

Jedyny selector-safe minimalny operator retardacji bez anchoru ma postac izotropowa:

`K_ret,iso^(1)(c;L) = cos(phi_ret(L;c)) * I_2`

Wnioski:
- `c` kontroluje wspolna faze/opoznienie,
- ale operator pozostaje proporcjonalny do `I_2`,
- wiec nie rozbija residualnego `O(2)`.

Czyli bez `psi0` minimalny jawny operator `c`-based jest dla selektora trywialny.

## Case B: psi0 imported

Jesli importujemy anchor `psi0`, dopuszczalny minimalny operator retardacji ma postac:

`K_ret,psi0^(1)(c;L_parallel,L_perp) = R(psi0) diag(cos(omega*L_parallel/c), cos(omega*L_perp/c)) R(-psi0)`

Wtedy:
- gdy `L_parallel = L_perp`, znow dostajemy operator izotropowy,
- gdy `L_parallel != L_perp`, pojawia sie rzeczywisty split widmowy/odpowiedzi.

Naturalna miara nowego efektu:

`Delta_ret = cos(omega*L_parallel/c) - cos(omega*L_perp/c)`

Jesli `Delta_ret != 0`, lane daje nowy pair-level efekt.

## Result

Najuczciwszy wniosek jest taki:

- samo `c` nie daje selector source,
- `c` moze sterowac szybkoscia/faza sprzezenia,
- nowy efekt pair-level pojawia sie dopiero po imporcie anisotropowej geometrii drogi albo anchoru `psi0`,
- wiec `c` jest sensownym skladnikiem operatora, ale nie samodzielnym rozwiazaniem problemu selektora.

## Frontier

`H42_B1 := a minimal c-based retardation operator on pair1 is selector-trivial without an imported anchor and can yield a genuine new spectral/response split only through psi0-dependent anisotropic path data, so c controls coupling speed but does not by itself generate selector orientation`

## Hard limits

- no `theorem-level PASS`
- no `full-closure PASS`
- no claim that `c` alone selects orientation
- no claim that `psi0 + c` is strict core
- no claim that `QW-2191` is discharged
