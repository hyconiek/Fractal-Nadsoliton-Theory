# Reproduction

Run `python verify.py` from this directory (SymPy and NumPy). The script reads no network or git resources, independently rebuilds two comparison jump matrices and writes RESULTS.json. Exact checks cover trace powers 1..12, the 4-partition identities and TV, Newton polynomial on rational and exact-radical data and a shell-ordering countermodel. A second floating enumeration iterates over 12^4 ordered quadruples with compensated summation. The code does not establish any physical model, independent source law or formal theorem certification.

Run `sha256sum -c MANIFEST.sha256` from this directory for integrity. Input hashes in INPUTS.sha256 name exact mounted source files that were available during execution; they are not a live repository checkout attestation.
