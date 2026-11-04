# Replay

From the repository root:

```sh
python3 fin_continuum_normal_form_research_20261104/verify_research.py
```

Python standard library only. The verifier writes verification.json.
Exact rational matrix controls check conjugation, idempotence, force
projection and full-generator commutation, with negative controls for
componentwise projection. The matrix example uses a real invertible
phase encoding instead of the complex half-wave encoding; conjugation
identities are the same, and it is not a spectral discretization of the PDE.
Other checks cover finite source evaluations, scaling identities, floating
homological-symbol previews and pinned input hashes. They do not formally
prove the infinite-dimensional estimates or simulate the actual trajectory.
Historical source/checker outputs remain unchanged.
