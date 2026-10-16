# Replay and source provenance

1. `python3 verify_package.py` validates the standalone package manifest and completeness.
2. `python3 verify_fp.py` runs 27 bounded checks and writes `FINITE_DIAGNOSTICS.json`. Its exact rational toy and high-precision non-certified integrals are **not** proof of q12 all-N feasibility.
3. `PROVENANCE.json` and `INPUTS.sha256` bind only actually accessible source files. Latest review seed and historical/current wrappers listed as unavailable are **not** silently substituted by predecessor archives.
4. Review exact all-N primal/dual statement in `FP-002/CONDUCTANCE_CONSTRUCTION_OR_DUAL.md`, necessary conditional cap theorem in `FP-002/POSITIVITY_AND_CAP_PROOF.md`, and open atom. There are no large count-lattice enumerations or external experiments.
