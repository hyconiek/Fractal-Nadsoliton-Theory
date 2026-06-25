# P3075/S2025 variational-source obstruction table

Status: `P3075_LOCAL_DIRICHLET_VARIATIONAL_SOURCE_PARTIAL_EXPORT_MEAN_CENTERING_NONLOCAL_OBSTRUCTION`

## Finite certificate
- content grep lanes: `4`
- content grep hits: `33888`
- P3074 accepted internal Lyapunov rows: `1008`
- accepted P3073 flow rows reused: `192`
- candidate generators: `5`
- variational-source matrix rows: `960`
- exact negative-gradient match rows: `192`
- accepted local variational-source rows: `96`
- local Dirichlet accepted rows: `96`
- global variance exact-but-nonlocal rows: `96`
- mean-centering local source rows: `0`
- satisfied proof obligations: `4/6`

## Decision
P3075 constructs the missing variational-source interface for the P3074 monotone flows.  The local Dirichlet quadratic generator has exact negative-gradient rows matching all 96 accepted cycle-Laplacian target rows, so a scoped internal local gradient source is exported for that flow.  The 96 mean-centering target rows match exactly only through the global variance generator, which is intrinsic but nonlocal on Z12, so no local source for mean-centering is exported.  The result remains finite, dimensionless, and internal: it is not a unit-bearing action, spacetime EOM, gauge field, or empirical physics map.

## Recommendation
Use only the P3075 local Dirichlet/Laplacian source as an internal dimensionless gradient-flow object.  The next proof-grade move is one bounded continuum-limit/spectral-dispersion audit: diagonalize the Z12 Dirichlet Laplacian source, compute its exact mode spectrum and small-k dispersion proxy, and test whether it has a lightlike/wave-compatible branch or only a diffusive/internal smoothing branch.  Do not promote the branch to observed light, gauge photons, spacetime EOM, units, or empirical physics without a separate unit/coordinate/source theorem.
