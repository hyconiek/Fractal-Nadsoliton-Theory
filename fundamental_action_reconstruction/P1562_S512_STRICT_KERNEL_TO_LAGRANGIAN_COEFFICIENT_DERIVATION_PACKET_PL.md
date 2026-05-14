# P1562 S512 Strict-Kernel To Lagrangian Coefficient Derivation Packet (No Legacy Bridge)

Status: `P1562_PROPOSED_STRICT_KERNEL_TO_LAGRANGIAN_COEFFICIENT_DERIVATION_PACKET`
As of: `2026-05-14`

## Cel

Wykonać realną pracę teoretyczno-obliczeniową bezpośrednio z kernela strict:

\[
K_{strict}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
\]

oraz wyprowadzić kandydaty współczynników do pełnego lagranżianu strict-only.

## Decyzja profesorska

Definiujemy momenty jądra:

\[
M_n=\int_0^{D} d^n K_{strict}(d)\,\mathrm{d}d,\quad n=0,1,2,3.
\]

Z nich budujemy bezwymiarowe wskaźniki:

\[
R_1=\frac{M_1}{M_0},\quad R_2=\frac{M_2}{M_0},\quad R_3=\frac{M_3}{M_0}.
\]

Następnie mapowanie do współczynników efektywnych (kandydaty strict-side):

\[
\lambda_{SM}^{(eff)} := |R_1|,
\quad
\kappa_{GR}^{(eff)} := |R_2-R_1^2|,
\quad
\epsilon_{mix}^{(eff)} := \frac{|R_3|}{1+|R_2|}.
\]

## Znaczenie fizyczne

- \(\lambda_{SM}^{(eff)}\): skala sprzężeń sektora materii,
- \(\kappa_{GR}^{(eff)}\): skala krzywiznowa sektora grawitacji,
- \(\epsilon_{mix}^{(eff)}\): skala mieszania SM–GR.

To nie jest import legacy — to strict-side derivation z `K_strict`.

## PASS/FAIL

PASS jeśli momenty i współczynniki są policzone stabilnie numerycznie
na zadanym oknie \([0,D]\) i eksportowane do artefaktu.

FAIL jeśli integracja niestabilna lub wartości niefizyczne (NaN/inf).

## Omówienie dla laika

To jak policzenie „odcisku palca” kernela strict:
zamiast deklarować współczynniki z góry, wyciągamy je bezpośrednio z matematyki
jądra i dopiero potem wkładamy do lagranżianu.
