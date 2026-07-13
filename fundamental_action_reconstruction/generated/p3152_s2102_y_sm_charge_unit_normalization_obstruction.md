# P3152/S2102 Y_SM charge-unit normalization obstruction

Status: `P3152_YSM_CHARGE_UNIT_NORMALIZATION_SCALE_TORSOR_NO_STRICT_SOURCE`

## Constructed object
- `U_Y^tors charge-unit normalization torsor`
- Classification: `bounded_scale_torsor_obstruction`
- Scope: `nonzero rational rescalings of the P3150 hypercharge ray, audited against Yukawa/anomaly residuals and local electric-charge values`

## Finite theorem
`P3152_T1_hypercharge_charge_unit_scale_torsor_obstruction`: The P3150 hypercharge equations and redundant anomaly checks are homogeneous in Y.  On the sampled nonzero rational scale orbit, every row has zero Yukawa/anomaly residuals while Y(H), hypercharge gcd, and electric-charge values vary.  Therefore current algebraic compatibility fixes only a projective hypercharge ray, not the physical charge unit or the convention Y(H)=1/2.

## Finite counts
- `sampled_nonzero_rational_scales`: `8`
- `invariant_scales`: `8`
- `distinct_Y_H_values`: `8`
- `candidate_source_rows`: `4`
- `strict_accepted_source_rows`: `0`

## Scale orbit rows
- `c=-3`: `Y(H)=-3/2`, residuals zero `True`, hypercharge gcd `1/2`, electric-charge gcd `1`
- `c=-2`: `Y(H)=-1`, residuals zero `True`, hypercharge gcd `1/3`, electric-charge gcd `1/6`
- `c=-1`: `Y(H)=-1/2`, residuals zero `True`, hypercharge gcd `1/6`, electric-charge gcd `1/3`
- `c=-1/2`: `Y(H)=-1/4`, residuals zero `True`, hypercharge gcd `1/12`, electric-charge gcd `1/12`
- `c=1/2`: `Y(H)=1/4`, residuals zero `True`, hypercharge gcd `1/12`, electric-charge gcd `1/12`
- `c=1`: `Y(H)=1/2`, residuals zero `True`, hypercharge gcd `1/6`, electric-charge gcd `1/3`
- `c=2`: `Y(H)=1`, residuals zero `True`, hypercharge gcd `1/3`, electric-charge gcd `1/6`
- `c=3`: `Y(H)=3/2`, residuals zero `True`, hypercharge gcd `1/2`, electric-charge gcd `1`

## Candidate source rows
- `Y_SM^ray homogeneous Yukawa/anomaly system`: fixes scale `False`, strict source `False`, noncircular `True`; equations select only a projective ray; every sampled nonzero rational scale remains valid
- `P3148/P3149 installed SM electric-charge convention`: fixes scale `True`, strict source `False`, noncircular `False`; sets Y(H)=1/2 by convention/installed registry, circular for strict source selection
- `P3146 axiom unit/action bridge`: fixes scale `False`, strict source `False`, noncircular `True`; length/time/action postulates do not define the U(1) charge unit
- `P3116 dimension-source functor inventory`: fixes scale `False`, strict source `False`, noncircular `True`; current physical-unit functors remain unsourced and do not export a charge valuation

## Decision
P3152 constructs the missing charge-unit object U_Y^tors and proves a finite scale-torsor obstruction for strict normalization.

## Why this is not strict
The only current row that fixes Y(H)=1/2 imports the installed SM charge convention.  Existing unit/action and dimension-source audits do not export a noncircular U(1) charge valuation from strict nadsoliton data.

## Recommendation
Pivot away from SM charge normalization unless a new strict charge-valuation object is introduced.  The next proof-grade route is P3153: a GR/EH nonproxy coupling audit for the axiom branch, testing one explicit metric/EH source interface without claiming L_total, SM generation, selector closure, bridge completion, role transfer, or ToE.
