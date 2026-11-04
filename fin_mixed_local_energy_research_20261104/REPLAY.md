# Replay

From the repository root:

```sh
python3 fin_mixed_local_energy_research_20261104/verify_research.py
```

Python standard library only. The checker separates exact polynomial and
rational controls, numerical integration-by-parts previews, and pinned
source hashes. It writes verification.json. The shifted-symbol polynomial
identity and Gaussian-kernel checks exercise the weight-transfer algebra;
they are not a discretization or formal verification of the PDE theorem.
The analytical proof supplies the uniform operator estimates. No historical
input or checker output is modified.
