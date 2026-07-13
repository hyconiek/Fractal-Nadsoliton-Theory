# P3156/S2106 alpha_geo to Higgs-coupling formula audit

Status: `P3156_ALPHA_GEO_TO_HIGGS_COUPLING_FORMULA_AUDIT_NO_STRICT_UNIT_SOURCE`

## Constructed object
- `Phi_alpha^H alpha_geo-to-Higgs-coupling formula schema`
- Classification: `single_formula_schema_dimension_and_source_audit`
- Scope: `tests whether S=alpha_geo can source both mu2 and lambda with units`

## Finite theorem
`P3156_T1_alpha_geo_higgs_coupling_unit_obstruction`: The single audited schema Phi_alpha^H can make dimensionless Higgs-coupling ratios from S=alpha_geo, but the Higgs mass parameter mu2 has mass dimension two.  Raw alpha_geo formulas are dimensionally invalid for mu2; dimensionally valid variants require an imported mass^2 scale M2.  Current artifacts therefore do not export a strict scalar-to-Higgs-coupling source for both mu2 and lambda with units.

## Finite counts
- `formula_schema_count`: `1`
- `formula_variant_rows`: `3`
- `dimensionally_valid_rows`: `2`
- `rows_importing_mass_scale`: `2`
- `accepted_strict_formula_rows`: `0`

## Formula rows
- `raw_dimensionless_alpha`: mu2 `S` dim `0`, lambda `S` dim `0`, v2 `1`, valid `False`, imports M2 `False`, accepted `False`
- `alpha_ratio_with_imported_mass_scale`: mu2 `M2*S` dim `2`, lambda `S` dim `0`, v2 `M2`, valid `True`, imports M2 `True`, accepted `False`
- `alpha_squared_lambda_with_imported_mass_scale`: mu2 `M2*S` dim `2`, lambda `S**2` dim `0`, v2 `M2/S`, valid `True`, imports M2 `True`, accepted `False`

## Decision
P3156 satisfies the P3155 continuation rule by auditing exactly one new formula schema.  It fails as a strict source because alpha_geo is dimensionless and the only dimensionally valid variants import M2.

## Recommendation
Do not continue the axiom-branch Higgs/SM/EH lane without a genuinely new strict mass-unit source.  The next proof-grade move should pivot to the independent Omega_dim/K_dim unit-source obligation from P3116, or to a selector-source intake if a new non-premise source object is supplied.
