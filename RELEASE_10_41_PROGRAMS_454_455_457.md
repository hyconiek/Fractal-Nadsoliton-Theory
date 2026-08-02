# Release 10.41 — FIN Local Research Programs P454–P457

## A Certified Nested Comb Dual, an Ordered-Symmetry Obstruction, and a Six-Decimal Global Cover

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-08-01  
**Language:** English  
**License:** CC BY 4.0

## Abstract

Release 10.41 records three local-only analytical and computer-assisted FIN
programs. P454 derives the exact nested semidefinite dual of the declared
three-slot causal-discrimination problem and supplies an explicit rational
dual certificate. Together with the inherited P451 primal witness, it bounds
the global full-cone optimum to a gap below `6.6e-6`. P455 proves that realness
and global bit-complement symmetry leave a five-dimensional fixed space, two
dimensions larger than the O167 coherent ansatz, thereby refuting a proposed
symmetry-only reduction. P457 uses the P452 symmetrization theorem and a new
outward interval cover to confine the declared full-simplex coarse-erasure
optimum to a gap below `1e-6`.

## Principal results

- Global three-slot causal half-distance bracket:
  `[0.52332810026048937, 0.523334700252]`.
- Exact bracket width:
  `0.00000659999151063`.
- Corresponding success-probability bracket:
  `[0.761664050130244685, 0.761667350126]`.
- Exact nested dual ladder of sizes `8 -> 4 -> 2 -> 1`, with all six matrix
  slacks and both scalar slacks certified positive.
- Certified minimum trigonometric dual-slack lower bound:
  `4.9999999969485233e-7`.
- Exact real complement-fixed causal dimension: `5`; exact O167 face
  dimension: `3`; symmetry-allowed residual dimension: `2`.
- Eight deterministic five-dimensional searches produce no improvement over
  the O167 value above `6.7e-16`; this remains strong evidence, not a theorem.
- Full-simplex coarse-erasure bracket:
  `[0.463278283192093, 0.46327928294340853]`.
- Coarse-erasure gap: `9.997513155113325e-7`, approximately `992.7` times
  narrower than the Release-10.40 bound.

## New typed objects

- O170 — Nested Comb Dual Ladder.
- O171 — Ordered-Comb Symmetry Residual Space.
- O172 — Six-Decimal Coarse-Erasure Cover.

## Scope and nonclaims

All work is local and dimensionless. No laboratory record, external audit,
internet result, remote computation, or physical evidence is used. The exact
P454 optimizer and its uniqueness remain open. P455 proves that the known
symmetries do not force O167; its numerical non-improvement result is not an
optimality theorem. P457 retains the exact P452 code, parameter, and heralded-
loss scope. The legacy and strict kernels remain separated; no legacy
physical role is transferred. QW-2191, canonical orientation, dimensional
units, apparatus, complete legacy-to-strict bridge, Standard Model, gravity,
`L_total`, and ToE closure remain open or unclaimed.

## Main document

`FIN_Local_Research_Checkpoint_P454_P457_EN.pdf`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1038 python3 fin_programs_454_455_457.py
python3 fin_programs_454_455_457_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P454_P457_EN.tex
python3 -m unittest -v test_fin_programs_454_455_457.py test_fin_checkpoint_p454_p457.py
sha256sum -c FIN_PROGRAMS_454_455_457_RELEASE_MANIFEST.sha256
```

## Next selected batch

P464 — close or sharply reduce the remaining P454 primal–dual gap through a
restricted dual/KKT inclusion; P458 — interval derivative isolation and
uniqueness test for P457; P459 — propagate P453 uniqueness through the O163
detector uncertainty box.

