# Replay

From the repository root with Python and SymPy (recorded: 1.14.0):

    python fin_virial_absorption_audit_20261107/check_virial.py

The script checks pinned inputs and writes verification.json. Identities
use exact symbolic simplification, not floating-point tolerances. The
49 checks comprise 31 source hashes and 18 analytical bookkeeping controls.

MANIFEST.json covers all other files in this package. Historical sources
are read only. PASS does not certify the analytic function-space extensions
or independently validate the complete conditional action theorem.
