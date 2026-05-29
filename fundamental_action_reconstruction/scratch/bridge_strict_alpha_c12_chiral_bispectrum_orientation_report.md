# Scratch strict-alpha C12 chiral bispectrum orientation probe

Status: C12-invariant / reflection-chiral Fourier phase audit; no strict selector discharge.

- Rows: `24`; separating bispectrum pairs: `['1,1', '1,2', '1,5', '2,3', '5,5']`.
- All selected pairs translation-invariant over sources: `True`.
- Source degeneracy per orientation: `{'-1': {'signature_count': 1, 'source_counts_per_signature': [12], 'representative_sources': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]}, '1': {'signature_count': 1, 'source_counts_per_signature': [12], 'representative_sources': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]}}`.
- Target replay: `q^5=256/243`, eta residual `0.000e+00`.
- Honest read: chiral bispectrum phase can distinguish orientation only after handedness is allowed, and it leaves source unresolved.
- No false pass: no strict chirality/source theorem, no QW-2191 discharge, no ToE closure.
