# P2687/S1637 one-new strict anisotropic source-class audit

Status: `P2687_ONE_NEW_STRICT_ANISOTROPIC_SOURCE_CLASS_AUDIT_NO_FALSE_PASS`

## Content-first grep
- `p2686_continuation_rule`: `56` hits
- `derived_lapse_energy_source`: `75` hits
- `non_energy_neutral_tensor_transport`: `119` hits
- `provider_nonexport_registry`: `15` hits
- `forbidden_imports`: `16858` hits

## Source-class candidate rows
- `derived_lapse_energy_source` targets `P1977_positive_energy_no_go`: exported=`False`, evades_now=`False`, verdict=`NEAR_MISS_NOT_EXPORTED`; requirement=`rho_base + u0 split must satisfy u0 = sigma1**2 + sigma1*sigma2 + sigma2**2 to offset rho_required = -sigma1**2 - sigma1*sigma2 - sigma2**2`
- `non_energy_neutral_tensorial_shear_transport` targets `P1978_energy_neutral_transport_obstruction`: exported=`False`, evades_now=`False`, verdict=`NEAR_MISS_NOT_EXPORTED`; requirement=`transport must abandon U00=0 and carry required_U00 = -sigma1**2 - sigma1*sigma2 - sigma2**2`

## Obligation matrix
- `derived_lapse_energy_source` / `typed strict source object`: `False`
- `derived_lapse_energy_source` / `internal derivation/provider citation`: `False`
- `derived_lapse_energy_source` / `evade P1977_positive_energy_no_go`: `False`
- `non_energy_neutral_tensorial_shear_transport` / `typed strict source object`: `False`
- `non_energy_neutral_tensorial_shear_transport` / `internal derivation/provider citation`: `False`
- `non_energy_neutral_tensorial_shear_transport` / `evade P1978_energy_neutral_transport_obstruction`: `False`

## Verdict
P2687 audits the exact P2686 continuation: one new strict anisotropic source class.  The two finite candidates are explicit and computable, but neither is currently exported as a strict typed source/provider.  A lapse/energy split would need an internal source with u0=Q_shear to avoid the P1977 negative-rho conclusion; a non-energy-neutral tensorial transport would need to carry the required U00 component to evade P1978.  Both are near-miss design equations, not current theorem exports, so this lane freezes as bounded no-go rather than reopening reverse closure.

## Next honest step
P2688 should pivot back to the broad state-map and choose a different live proof frontier, with first preference returning to the kernel-split-robust P46/N49 direct-route zero-witness/no-go matrix named by F3/P2681 (especially the direct m2 psi4 target-action coefficient defect on common psi4**2/2 support).  Do not replay Lagrangian/EOM reverse closure unless a genuinely exported strict anisotropic source class is introduced.
