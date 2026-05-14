# P1563 S513 Strict Kernel -> Coefficients -> Lagrangian -> Euler-Lagrange EOM Export Packet (No Legacy Bridge)

Status: `P1563_PROPOSED_STRICT_KERNEL_TO_EOM_EXPORT_PACKET`
As of: `2026-05-14`

## Cel

Domknąć pełny tor wymagany przez Ciebie:

`kernel strict -> współczynniki -> lagranżian -> równania ruchu`.

## Wejście (z P1562)

Z `K_strict(d)` otrzymane współczynniki:

- `lambda_sm_eff`,
- `kappa_gr_eff`,
- `epsilon_mix_eff`.

## Konstrukcja lagranżianu strict-side

Przyjmujemy reprezentację roboczą:

\[
\mathcal{L}_{total}=
\underbrace{\tfrac12(\partial\psi)^2-\tfrac12\lambda_{SM}^{(eff)}\psi^2}_{\mathcal{L}_{SM}}
+
\underbrace{\kappa_{GR}^{(eff)}R}_{\mathcal{L}_{GR}}
+
\underbrace{\epsilon_{mix}^{(eff)}\psi R}_{\mathcal{L}_{mix}}.
\]

## Równania ruchu (Euler–Lagrange)

Dla pola \(\psi\):

\[
\Box\psi + \lambda_{SM}^{(eff)}\psi - \epsilon_{mix}^{(eff)}R = 0.
\]

Dla sektora geometrycznego (schematycznie):

\[
\kappa_{GR}^{(eff)} G_{\mu\nu}
=
T^{(\psi)}_{\mu\nu}
+
\epsilon_{mix}^{(eff)}\Theta^{(mix)}_{\mu\nu}.
\]

To jest strict-only eksport EOM bez legacy bridge.

## PASS/FAIL

PASS = formalny eksport \(\mathcal{L}_{total}\) i obu EOM do artefaktu.

FAIL = brak któregokolwiek członu toru kernel->coeff->L->EOM.

## Omówienie dla laika

To już kompletna „linia produkcyjna”:
z kernela liczymy liczby, z liczb budujemy wzór energii układu,
a z tego wzoru wyprowadzamy równania, które mówią jak układ się zachowuje.
