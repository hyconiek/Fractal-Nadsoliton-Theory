# SRN-only replay

From this directory run `python verify_srn.py` to regenerate `INDEPENDENT_VERIFICATION.json` and check exact rational projector identities, finite generator jumps, conditional hypergeometric fixtures, raw gauge invariance, nominal ideal angular calculations and exact-root exponent bookkeeping. Run `python finalize_srn.py` to regenerate input hashes, a comprehensive SHA-256 manifest and the SRN-only ZIP. Inputs are hashed from actual available bytes and not copied into the package.

Requirements: Python 3, numpy, sympy, mpmath. Tests are algebra/numeric diagnostics, not mechanized proofs of local Gibbs/MT or asymptotic all-N results. Read each stage's analytical proof separately. The source user plan and external archives are input files, not executed commands. Absence of a source review file must not be represented as a verified historical artifact.
