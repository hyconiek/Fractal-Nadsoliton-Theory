# Replay

From the repository root:

```sh
python3 fin_local_feedback_research_20261104/verify_research.py
```

Standard-library Python. The verifier checks the exact polynomial split,
finite rational absorption/scaling controls, numerical weighted-integration
identities and pinned source hashes. It writes verification.json. None of
these checks is a PDE simulation, formal proof or independent review.
The packet countercontrol is analytical; the finite support checks do not
simulate the actual prepared trajectory. The polynomial helper code is
reused from QH without executing or modifying its checker or outputs.
