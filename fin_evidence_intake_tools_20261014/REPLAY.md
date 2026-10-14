# Replay

From this directory, with Python, numpy, scipy and mpmath available:

```bash
python run_unit_tests.py
python independent_verify.py
python strict_records.py fixtures/valid.csv
python strict_records.py fixtures/nan_elapsed.csv  # expected exit 2
python strict_process_gate.py --model FIN_A7 --kinetic heat_bath --gain 3 --shortcut
```

For a future decision package, provide a JSON object `{"n": <positive integer>, "counts": [13 exact integer counts]}` to `strict_decision.py`; use only after the calibration and sampling premises have independently passed.

Do not treat these commands as a physical experiment.
