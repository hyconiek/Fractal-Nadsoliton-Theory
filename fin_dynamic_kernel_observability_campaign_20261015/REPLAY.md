# DKO bounded independent replay and manifest verification

1. With Python 3, numpy and mpmath installed, run `python3 independent_verify.py` inside this new DKO directory. Script rebuilds the source kernel from declared omega/phi/eta, exact rational P3/P4, the actual q12 per-slot LOO-HB rates, the finite generator brackets for N=120, and limiting Gaussian sign and Toeplitz KL fixtures. It creates `INDEPENDENT_VERIFICATION.json`.
2. Run `python3 verify_package.py`; it reads MANIFEST.sha256 and checks every manifest entry and stage required files. The manifest excludes its own bytes and PACKAGE_VERIFICATION.json to avoid circular hashing.
3. The analytic proof is contained in DKO-001/ACTUAL_CHARACTERISTICS.md, DKO-002/JOINT_FAST_PATH_THEOREM_OR_BLOCKER.md, DKO-003/BOUNDED_SIGN_LIMIT.md and CLOCK_AND_AMPLITUDE_GAUGE_THEOREM.md, and DKO-004/DEPENDENT_RECORD_INFORMATION.md / LOWER_AND_UPPER_DECISION_BOUNDS.md. Re-running finite assertions does NOT mechanically validate the infinite-N martingale-problem argument.
4. SOURCE_INPUTS.sha256 hashes the real available supplied plan and predecessor archive bytes. The KO intake verdict pasted in chat and inaccessible source review files are deliberately not given invented hash values.
5. No fresh physical data acquisition, web source correction, repository modification or empirical test occurs in this replay.
