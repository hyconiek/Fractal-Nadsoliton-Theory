# D2-SOURCE-01/02/03 — marginal point and no-go

Status: **NO_GO_PLUS_CONDITIONAL_SELECTOR**

## Main result
Locality, positivity, exact refinement and power-law self-similarity do not select D_H=2: every r>1 has a valid local bulk.  For regular b-ary trees the local conductance runs as `(r^2/b)^d`.  Only the stronger condition of unrescaled uniform nondegeneracy (or zero scaling dimension of the local conductance) forces `r^2=b`, hence `D_H=2`.  The repository does not currently source that stronger condition.

## Core formulas
\[
g_d=rac{r^2}{b(r^2-1)}(r^2/b)^d.
\]
Uniform `0<g_-<=g_d<=g_+<infinity` for all d implies \[
r^2=b,\qquad D_H=2.
\]

## Evidence / reproduction
Repository searches found self-similarity/locality premises explicitly marked unsourced; ST753 also leaves scale-charged data as the escape hatch.

## Caveats
`D_H=2` is not a strict FIN dimension result.  Natural storage `M_d=b^-d` already gives `tau_d~ell_d^2` for every D_H.

## Next question
Keep D2-SOURCE stopped unless a genuinely new scale-charged strict datum appears.
