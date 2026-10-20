# Reproduction

With numpy/scipy/sympy available, execute `python3 verify.py` in this directory. The script regenerates RESULTS.json and runs exact symbolic, Markov spectral, counterfamily and finite-noise checks. Then run `sha256sum -c MANIFEST.sha256`.

`INPUTS.sha256` records the exact already-uploaded source snapshots at paths relative to the research directory inside `/mnt/data`. Input verification (`sha256sum -c INPUTS.sha256`) requires keeping the source uploads at those relative locations. An extracted ZIP alone has no authority to claim access to the original uploaded AGENTS or full live repository.

No stochastic path sampling, calibration, source selection, or peer-reviewed formalization is performed.
