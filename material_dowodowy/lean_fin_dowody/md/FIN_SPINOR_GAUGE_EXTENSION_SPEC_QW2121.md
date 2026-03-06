# FIN Spinor+Gauge Extension Spec (QW-2121)

This document is a formal extension specification.
It does not claim strict-derivation status for spinor/gauge/gravity action blocks.

## Core Extended Action (symbolic)

L_total = L_scalar + L_psi + L_Y + L_gauge + L_gravity_bridge

### L_scalar (existing strict chain block)
- L_scalar = sum_o [1/2 d_mu Psi_o^dag d^mu Psi_o - V_o(|Psi_o|^2)] - 1/2 sum_{o!=o'} K_{oo'} Psi_o^dag Psi_{o'}

### L_psi (spinor kinetic extension)
- L_psi = sum_i bar(psi_i) (i gamma^mu D_mu - m_i) psi_i

### L_Y (Yukawa extension)
- L_Y = - sum_{ij} [ y^u_ij bar(Q_i) tilde(H) u_j + y^d_ij bar(Q_i) H d_j + y^e_ij bar(L_i) H e_j ] + h.c.

### L_gauge + D_mu (gauge extension)
- L_gauge = -1/4 G^a_{mu nu} G^{a mu nu} -1/4 W^i_{mu nu} W^{i mu nu} -1/4 B_{mu nu} B^{mu nu}
- D_mu = partial_mu - i g_s T^a G^a_mu - i g tau^i W^i_mu - i g' Y B_mu

### Gravity bridge statement
- G_{mu nu}(g_eff) = 8 pi G_eff T_{mu nu}^{eff} + Delta_{mu nu}^{UV}

## Status tags
- scalar block: strict-chain existing
- spinor/gauge blocks: formal extension, derivation pending
- gravity action-level derivation: pending
