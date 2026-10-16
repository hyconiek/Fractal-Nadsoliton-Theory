# HM-001: Exact cavity update

**Scientific disposition:** `EXACT_FINITE_N_IDENTITY`.

Completing the square gives the exact native LOO softmax after eliminating a freshly sampled cavity field. θ|n and ν_r are different; no (N−1) denominator or allowed-label renormalization.

See the dedicated deliverables in this task, the complete derivations in [HM_THEOREM_PROOFS.md](../HM_THEOREM_PROOFS.md), and the source-access disclosures in [INPUT_PROVENANCE.json](../INPUT_PROVENANCE.json). The controls in `verify_hm.py` are bounded convention checks rather than a mechanical proof of the growing-N limits.
