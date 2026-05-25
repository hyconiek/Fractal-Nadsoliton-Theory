# P2028 S978 Strict B1 GB-Quotient Counterterm Identifiability Theorem

Status: `OPEN_OBSTRUCTION_WITH_TRACE`

Local verdict: `PASS_QUOTIENT_THEOREM_WITH_TRACE`

P2027 showed that the strict B1 scalar operator profiles have rank `3` with
one Gauss-Bonnet null direction:

`n_GB = (1,-4,1,-1)`.

P2028 exports the quotient map

`T(a_R2,a_Ric2,a_Riem2,a_GB) = (a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB)`.

The canonical section is

`s(q_R2,q_Ric2,q_Riem2) = (q_R2,q_Ric2,q_Riem2,0)`.

Exact checks:
- `T(n_GB)=0`
- `T*s=I_3`
- `s*T` is idempotent
- `s*T(n_GB)=0`

Target quotient coefficients:
- R2_bar=9.9999999999999922e-01
- Ric2_bar=1.1656308464946203e-15
- Riem2_bar=6.0663244882429037e-17

The whole family

`a(tau)=s(T(a)) + tau*(1,-4,1,-1)`

has the same B1 scalar quotient class.  Across sampled tau values the maximum
full-system residual L-infinity norm is `4.219e-15`, with
tolerance `1.000e-10`.

## Honest Interpretation

This is a local quotient theorem.  It licenses the rank-3 B1 scalar
counterterm class modulo the GB null direction.  It does not license an
independent `a_GB`, four-channel uniqueness, a tensor-resolved projection,
background-global renormalization, unitarity, selector closure, or ToE closure.
