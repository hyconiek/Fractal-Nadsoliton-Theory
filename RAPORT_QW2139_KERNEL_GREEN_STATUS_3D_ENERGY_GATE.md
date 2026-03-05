# RAPORT QW-2139: KERNEL GREEN-STATUS + 3D ENERGY GATE

- Date UTC: 2026-03-04T19:44:20.197163+00:00
- Verdict: **KERNEL_GREEN_STATUS_3D_ENERGY_GATE_PASS_PARTIAL_ROLE_CLARIFIED**
- pass_count: `9/10`

## Kernel
- formula: `K(d)=cos(omega*d+phi)/(1+beta*d^eta)`
- omega=0.18575, phi=0.1625, beta=1.0, eta=1.8

## Green-status diagnostics
- Laplace residual ratio: `6.740236`
- Helmholtz residual ratio (best lambda): `4.801738`
- best lambda: `4.730126`

## 3D integrability diagnostics
- abs integral slope (log-log): `1.198604` (expected `1.200000`)
- L2 tail last increment: `4.788202e-03`
- grad tail last increment: `1.758524e-04`

## Scope boundary
- Full constructive Green-operator derivation from FIN action: `False`
