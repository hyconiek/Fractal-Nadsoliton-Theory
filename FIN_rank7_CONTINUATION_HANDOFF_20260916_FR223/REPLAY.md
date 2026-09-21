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

## 2026-09-15 continuation frontier replay

Run from the bundle root with bytecode disabled:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src python verify.py
```

Then replay durable continuation tests as fresh shards:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR1_residual_tails.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR8_FR9_staircase.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR10_diagonal_uv.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR11_low_e_r_arm.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR12_pure_r_arm.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR13_diagonal_bulk.py
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src pytest -q tests/test_FR14_strengthened_rv.py
```

Expected continuation total: **21 passed, 0 failed**.

FR2/FR3 are preserved proof notes but do not yet have standalone replay scripts in this continuation tree. FR4--FR7 are explicitly marked replay-required in `CONTINUATION_LEDGER.json`.
