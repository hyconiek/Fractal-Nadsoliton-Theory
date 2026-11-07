# Replay

From the repository root with Python and SymPy (recorded version 1.14.0):

    python fin_action_law_audit_20261107/check_audit.py

The script checks 30 pinned source hashes and 11 exact controls, then
writes verification.json. The controls cover projection conjugation,
incorrect operator ordering, finite-primitive endpoints, carrier inventory
and critical-time exponents. No floating-point tolerance is used.

MANIFEST.json covers all other package files. Historical inputs are read
without modification. The result is a bounded algebra/integrity audit,
not an independently reviewed mathematical theorem.
