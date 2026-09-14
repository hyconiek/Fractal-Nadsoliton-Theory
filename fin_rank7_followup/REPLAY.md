# Portable replay instructions

Run from the `fin_rank7_followup` directory. Python dependencies are recorded in `environment.json`.

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=.:..:src
python verify.py
python src/final_audit.py
python -m pytest -q tests/test_schema_and_resume.py
python -m pytest -q tests/test_G_checker.py
python -m pytest -q tests/test_I_local_cone.py
```

The complete campaign regression is intentionally sharded by test file. The recorded final result is 94 passed, 0 failed, 0 skipped; see `verification_final.json`. A monolithic pytest process exceeded the execution cap and is not the canonical replay mode.

For full sharded replay:

```bash
python replay_sharded.py
```

The historical source ZIP listed in `NONREPLAYED_INPUTS.json` is not required to run the promoted campaign package; it is explicitly not claimed replayed. Do not add or regenerate it merely to silence provenance reporting.
