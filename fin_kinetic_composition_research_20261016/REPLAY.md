# Replay and provenance

From the repository root:

```bash
python3 fin_kinetic_composition_research_20261016/verify_research.py
sha256sum -c fin_kinetic_composition_research_20261016/INPUTS.sha256
```

From this directory:

```bash
sha256sum -c MANIFEST.sha256
```

`verify_research.py` uses only Python's standard library. All decisions use
exact Fractions. Float outputs are previews only. It regenerates RESULTS.json
deterministically and checks:

- rational q12 tilted-Hamiltonian symmetry and Hessian rank, with a nonzero-current negative control;
- seven-state symmetry, detailed balance, stationary marginal and exact generator intertwining;
- distinct relative rates at the same energy/equilibrium;
- the native rate-two eigenmode and a rigorous alternating-series enclosure for a regular conservative control.

These finite checks do not constitute a proof of all-g path convergence,
independent verification of the source literature, or empirical validation.
Analytical proofs and their assumptions are given separately.

Three review manifests were freshly checked: MM 14/14, HM 13/13, NT 13/13,
zero missing and zero mismatching entries. See
`evidence/review_manifest_check.json`. Historical campaign replays were not
rerun by this comparative study. Their previously admitted results remain
inherited, with their original limitations.

`INPUTS.sha256` binds actual repository files consulted, using paths relative
to the root. `MANIFEST.sha256` binds this submission's files except itself
and bytecode caches. Source captures are rendered browser text and extracted
sections, not original PDFs. Captured search pages are discovery records only.
The MPR text companion has trailing whitespace normalized for repository
hygiene; its original DOM text is preserved in `evidence/gradient_ldp_browser_raw.json`.

At the start, AGENTS and NT/HM/MM materials were already modified/untracked.
During this work an external commit `761e0731cb4e45669605611ef0c3b38b58a0eff8`
incorporated them and some intermediate literature captures. This agent did
not make that commit. It subsequently corrected its own MPR capture from a
nonexistent v3 page to the verified v1 text. Do not infer that all captures
were complete at the intermediate commit. Current bytes are bound here.

No historical source/review, AGENTS, kernel, simulator or running research
campaign was modified by this study. No external experiment, spending or
messaging was performed. Physical and scientific admission is independent
of syntax/hash PASS.
