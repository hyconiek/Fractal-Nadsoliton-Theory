# Active 30-round FIN research goal

New goal checkpoint, starting after ST8590. Reports are Markdown; no PDF
generation and no document-production pipeline. Previous completed studies
are background, not counted again toward these thirty rounds.

Authoritative baseline consulted: current AGENTS.md, the actual worktree,
and the replication-consistency proofs/code. The previous work changed
scientific evidence and is classified as progress, not a no-progress turn.

## Completed, with proofs in REPORT.md

1. **ST8591:** transition-resolved replication Gram identity; reconstructing
   G2 from Q and the full transition tensor. Next: test sufficiency.
2. **ST8592:** exact counterexample: a nonnegative Gram transition matrix,
   valid reversible strict two-copy generator and all graph symmetries can
   still have no three-copy extension. Next: determine whether supplying an
   actually extendible full G2 closes the hierarchy.
3. **ST8593:** two infinitely extendible families with exactly the same full
   G1/G2 but different G3. Next: replace the pair cutoff by any finite cutoff.
4. **ST8594:** exact rational moment construction matching every generator
   through any prescribed finite order and differing at the next order.
   Next: identify what all orders determine, and whether a finitely specified
   self-similar source law can close the infinite hierarchy without hiding
   an added physical postulate.

5. **ST8595:** all ideal population orders identify μ on a supplied compact
   noise curve; exact triangular moment reconstruction. Next: stability.
6. **ST8596:** the alias laws have TV=1 but exact W1=2/(m+1); an operational
   transition-TV bound prevents confusing nonidentification with uselessness.
7. **ST8597:** exact algebraic noise amplification; the result is not promoted
   to a minimax impossibility theorem outside feasible moment data.
8. **ST8598:** a specified binary affine recursion uniquely supplies an entire
   hierarchy, and its contraction is fixed by a supplied variance.
9. **ST8599:** a ternary innovation recursion shares that variance but changes
   fourth moments; unspecified self-similarity does not close the law.
10. **ST8600:** graph reflection acts trivially on this radial noise coordinate,
    so spatial symmetry does not source the fair-sign innovation assumption.

## Remaining requirements

Rounds 11–30 are not completed. Continue adaptively from the above evidence.
Maintain explicit legacy/strict and mathematical/physical boundaries.
Revisit relevant source classes across the repository when an actual new
typed candidate becomes available; do not repeat closed source inventories.
The final report must include all thirty completed rounds, proofs and
falsification outcomes, conditional premises, source coverage and the best
remaining route toward physical predictivity. No completion claim yet.

## Replay of this checkpoint

```sh
python3 fin_transition_hierarchy/research.py
python3 -m unittest discover -s fin_transition_hierarchy -p 'test_*.py' -v
```

The current suite has 18 scientific tests, including exact transition-count
Gram corrections for different targets from the same origin. No PDF,
archive or publication-layout work is part of this goal checkpoint.
