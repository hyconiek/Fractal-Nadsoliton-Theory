# Replay instructions

Python 3 with numpy/scipy available. Working directory: archive root `fin_critical_metastability_campaign_20261014`.

```bash
python mt_analysis.py
python independent_verify.py
sha256sum -c MANIFEST.sha256
```

`mt_analysis.py` imports only the versioned `source_context/theory_core.py` for a small *numerical* model check. `independent_verify.py` reconstructs the FIN Fourier matrix independently and checks exact rational boundary inequalities using stdlib Fraction. Both output JSON; exact same-environment replay is summarized in `CLEAN_REPLAY_VERIFICATION.json`. The source documents provide the analytical arguments separately. The manifest covers every file except itself; source archives and the two user-uploaded markdown sources are preserved verbatim.

No external experiment, global trajectory simulation, large type-lattice enumeration or automatic AGENTS edit is performed. Scientific acceptance is reserved for independent intake.
