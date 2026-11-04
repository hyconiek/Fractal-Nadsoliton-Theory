# Replay and scope

From the repository root run:

```sh
python3 fin_quadratic_modal_history_research_20261104/verify_research.py
```

Python standard library only. The output separates exact finite controls,
floating-point identity previews, and source hash checks. The script writes
verification.json; historical inputs remain unchanged. It does not perform
a PDE simulation or formally verify the analytical proof. INPUTS.json binds
the instruction context and analytical dependencies at the stated base.
