# UE replay / scope

Run on Python 3 with SymPy:

```bash
python3 REPLAY/verify.py
```

It re-derives 11 **exact arithmetic** checks for the 2×2 noncommuting negative-direction countercontrol, transfer-energy ratio and coarse-acceleration separation and same-field two-bump composition cross term. The script writes `REPLAY/exact_controls.json` and stops with an assertion failure if any identity changes. Fractions/integer matrices are exact; no digit-spectrum campaign, physical parameter fit or QW-program execution is performed. The mathematical proofs and their assumptions are in the named UE files, not replaced by these checks.

`INPUTS.sha256` identifies the one available materialized Library AGENTS snapshot and the campaign inputs created in this workspace. `MANIFEST.sha256` hashes the final generated output corpus; files that could not be downloaded are listed explicitly as `UNAVAILABLE` in `SOURCE_ACCESS.md` and are **never assigned invented hashes**.
