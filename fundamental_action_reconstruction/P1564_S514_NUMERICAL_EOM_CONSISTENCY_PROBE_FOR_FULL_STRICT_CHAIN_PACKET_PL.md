# P1564 S514 Numerical EOM Consistency Probe For Full Strict Chain Packet (No Legacy Bridge)

Status: `P1564_PROPOSED_NUMERICAL_EOM_CONSISTENCY_PACKET`
As of: `2026-05-14`

## Cel

Po symbolicznym eksporcie EOM (`P1563`) wykonać krok numeryczny:

- sprawdzić wewnętrzną zgodność równań ruchu,
- potwierdzić brak sprzeczności znaków i skal,
- utrzymać strict-only tor `F_Nadsoliton => L_SM + L_GR`.

## Metoda

Dla siatki próbnej \(d \in [0,D]\):

1. obliczamy \(R(d)\) z `K_strict` (proxy curvature),
2. podstawiamy do równania \(\psi\)-sektora,
3. liczymy residuum:

\[
\mathcal{E}_{\psi}(d)=\Box\psi(d)+\lambda_{SM}^{(eff)}\psi(d)-\epsilon_{mix}^{(eff)}R(d).
\]

PASS gdy norma \(\|\mathcal{E}_{\psi}\|_2\) jest poniżej progu.

## Omówienie dla laika

To test, czy nasze równania „trzymają się kupy” na liczbach,
a nie tylko na papierze.
