# MM-only reproducibility

At the directory of this MM campaign, run `PYTHONDONTWRITEBYTECODE=1 python3 verify_mm.py`, and check `verification.json`. The verifier uses `sympy`, `mpmath.iv`, and `numpy`; it performs exact polynomial/Lyapunov tests, fixed-state labelled jump sums, nominal spectral Gaussian calculations and outward arithmetic conditional on gamma3 lying in the explicitly declared interval. It does NOT simulate trajectories or prove the asymptotic martingale limit.

Run `python3 verify_package.py` to check manifest content coverage and hashes. No program rewrites AGENTS, the HM/MT/KO/DKO source archives, or repository state. `build_mm_documents.py` regenerates narrative scaffolding; its output may be reformatted but the scientific assumptions and checks are in the files. This file is not a script for reconstructing missing original uploaded wrappers or reviews.
