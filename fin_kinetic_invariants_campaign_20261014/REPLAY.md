# Replay — FIN kinetic invariants

Python 3 + numpy + scipy are required. All commands execute locally and use tiny 12D/11D/10D matrices; no large-N enumeration and no network:

```
python kinetic_analysis.py
python local_regularities.py
python independent_verify.py
python verify_package.py
```

The first producer reads a byte-preserved `source_context/theory_core.py` for numeric FIN reconstruction, evaluates the inherited fold nominally and computes the equilibrium three-rule mobilities and finite detailed-balance probes; the second probes a few positive/negative fold-line coordinates. The **independent** verifier uses a separate Fourier construction and a different fold equation, does not import producer/scientific reconstruction modules, and checks nonlinear Onsager identities.

Output may differ in the least significant bits across BLAS/scipy/numpy environments. `MANIFEST.sha256` binds the archived completed instance; running generators in-place may invalidate it. For unmodified-archive hash checks, run `sha256sum -c MANIFEST.sha256` before any replay, and make a copy for rerunning writers.

The mathematics in markdown remains an analytical proof with stated inherited premises; none of these commands proves a stochastic theorem or independently certifies mpmath.iv fold enclosures.
