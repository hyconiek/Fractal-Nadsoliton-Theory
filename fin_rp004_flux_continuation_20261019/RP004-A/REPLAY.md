# RP004-A replay

1. Run `python3 ../check_spectral_energy.py` from the parent campaign (or use its absolute path). This uses an independent two-root spectral energy quadrature and is **diagnostic only**.
2. Check `E/t` against `Omega*Gamma/2`; quadrature preview at a=-1/2: T40 0.12639336, T80 0.12825272, T160 0.12897544 versus limit 0.13028450.
3. Verify the exterior sign independently: `e_t+J_x=-S(1+cos Omega t)u_t`; on two exterior half-lines `F_u=E_ext(T)+int S(1+cos)u_t`.
4. Spectral/energy proofs are supplied in the markdown, not established by these bounded checks.
