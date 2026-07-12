# P3134/S2084 explicit legacy phase/torsion D_HL candidate audit

Status: `P3134_EXPLICIT_LEGACY_PHASE_TORSION_DHL_CANDIDATE_ORIGIN_POLARITY_NO_GO`

## Constructed object
`D_HL(r,k,lambda;x)=lambda*beta_tors*sin(2*pi*k*(x-r)/12)`

a beta_tors-scaled oriented phase-gradient/torsion residual centered at a support representative r with chiral polarity lambda

## Finite certificate
- beta_tors: `0.01`
- candidate defects: `120`
- nonzero odd candidates: `120`
- conditionally support-coupled candidates: `120`
- sampled translation symmetry rows: `240`
- rows equivariant when origin shifts: `240`
- quotient-invariant nonzero `t=1` rows: `0`
- accepted import-free D_HL sources: `0`

## Decision
P3134 constructs the missing object rather than only naming it. The explicit beta_tors-scaled sine residual gives 120 nonzero candidates that are odd around a chosen support representative and pair under lambda -> -lambda. This proves the legacy phase/torsion split can generate the right local mathematical shape for D_HL. The same computation also proves why it is not yet an import-free strict source: the formula is translation-covariant only if the origin r is carried along, not invariant after quotienting the diagonal Z12 origin orbit, and the polarity lambda remains paired. Thus the construction upgrades the gap from vague to precise: missing are an internal origin section and an internal lambda/sign law, not another generic helix label.

## Recommendation
Do not add another helical residual. The next proof-grade move is exactly one internal origin-and-polarity source test for the constructed D_HL: either export a strict law selecting (r,lambda), or prove a no-go for all nadsoliton-internal candidates currently available. A minimal computable target is a joint (r,lambda) selector-source matrix over the P3134 candidate family, testing only translation quotient, inversion polarity, Aut(Z12) behavior, and import freedom; no role transfer, bridge completion, L_total, or ToE promotion is licensed before that matrix passes.
