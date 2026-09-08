# Thirty-round FIN research campaign

This campaign started after ST8590. Reports are Markdown; no PDF
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
11. **ST8601:** identical recursive stationary law, different temporal
    generators, even after matching a clock using one observable.
12. **ST8602:** a valid hidden telegraph model has averaged drift Q but an
    extra C² in its second derivative, refuting homogeneous coarse closure.
13. **ST8603:** exact hidden-state elimination supplies a memory kernel and
    Schur response, conditional on the supplied environment rate.
14. **ST8604:** commuting fast-switching limit and a uniform operator bound;
    no finite-rate exact heat claim.
15. **ST8605:** extremal/flat moment data can identify special finite laws;
    this limits the interpretation of the generic finite-cutoff no-go.
16. **ST8606:** two-mode rigidity for positive common subordinators preserving
    the strict spectrum (up to a common calibration).
17. **ST8607:** a sharp clock-jump tail bound; tiny high-rate jumps remain
    indistinguishable from drift at finite precision.
18. **ST8608:** explicit positive one-mode counterexample tests necessity of
    two distinct nonzero spectral values.
19. **ST8609:** commuting stochastic maps can have negative operator
    eigenvalues and lie outside the positive heat-mixture class.
20. **ST8610:** full reading of K1/K2/F2/F3/S2/SUMMARY plus current guards;
    class-premise audit prevents falsely promoting clock rigidity to a
    sourced population composition or legacy/strict physical bridge.

21. **ST8611:** invisible identity attempts expose an observational gauge
    in common-event intensity measures. Next: identify the quotient law.
22. **ST8612:** all complete orders identify a finite intensity measure modulo
    its identity atom. Next: test whether the class contains independent dynamics.
23. **ST8613:** independent dynamics requires a separate component or a
    singular small-jump limit; a controlled tensor expansion verifies the limit.
24. **ST8614:** uniqueness of independent rates plus a common-event measure
    in the declared integrable class, with exact weighted-moment reconstruction.
25. **ST8615:** quantum common-phase channels match through any finite order
    and differ higher; unitary marginals provide an extremal exception.
26. **ST8616:** nondegenerate preparation/readout aliases survive after the
    dynamical law is fixed. Next: test kernel-split robustness.
27. **ST8617:** the signed-legacy cover supports exact pair-law aliases with
    different triple laws, without identifying legacy with strict.
28. **ST8618:** approximate self-similar source laws give rigorous predictive
    error bounds; finite bit truncation is checked exactly.
29. **ST8619:** internal-clock stopping gives a normalized Green operator and
    distinguishes a common scale convention from a changed prediction.
30. **ST8620:** adversarial synthesis, source coverage and completion audit;
    no claim of physical or ToE closure.

## Final audit

The thirty investigations and their proofs/results are documented in
REPORT.md. Completion of this campaign does not mean completion of FIN as
a physical theory. COMPLETION_AUDIT.md maps the actual research request to
evidence and states the remaining scientific limitations. Goal completion
requires that audit and the final test/replay checks to succeed.

## Replay of this checkpoint

```sh
python3 fin_transition_hierarchy/research.py
python3 fin_transition_hierarchy/temporal.py
python3 fin_transition_hierarchy/completion.py
python3 -m unittest discover -s fin_transition_hierarchy -p 'test_*.py' -v
python3 fin_transition_hierarchy/verify.py
```

The final suite has 46 scientific tests, including exact transition-count
Gram corrections for different targets from the same origin. No PDF,
archive or publication-layout work is part of this goal checkpoint.
