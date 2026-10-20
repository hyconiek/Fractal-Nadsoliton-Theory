# Replay instructions

Run locally with Python 3 and standard library only:

```
cd fin_complete_anonymous_history_campaign_20261020
python3 verify.py
python3 verify_package.py
```

`verify.py` performs exact rational `J`, `J²`, row stochasticity and the two-step lower bound with an independent short seeded first-appearance sample and evaluates the analytical L formula. It writes `RESULTS.json`. It does **not** verify the external Azuma-Hoeffding theorem or turn a finite sample into an almost-sure proof.

`verify_package.py` verifies hashes from `MANIFEST.sha256` and ensures `RESULTS.json` is present. Archive CRC check is separate and is done at packaging time. Full proof obligations are in `FULL_ANONYMOUS_IDENTIFIABILITY_OR_COUNTERPAIR.md` and `FINITE_PREFIX_RECONSTRUCTION_THEOREM.md`.

The user-provided original repo path is not mounted; do not represent this as an in-tree Git replay of the current checked-out branch.
