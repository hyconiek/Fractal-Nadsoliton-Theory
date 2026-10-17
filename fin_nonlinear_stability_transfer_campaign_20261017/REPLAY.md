# Replay and proof-check scope

Run `python3 REPLAY/verify.py` from the campaign directory (Python 3.13 with sympy). This recomputes `REPLAY/exact_checks.json` deterministically. Recorded here: 31/31 exact algebraic/small-matrix controls passed. It tests residual algebra for signed sextic, kink first integral, signed g plus/minus noncommuting control gaps, scalar phi4 comparison eigenfunctions, matrix commutators, direct-sum of 24 modes, modulation normalization and the sector work cancellation.

The proof of the all-time PDE orbital estimate rests on (i) minmax/Sturm/Loewner argument, (ii) implicit-function modulation, (iii) H1 local Taylor bound, (iv) conserved canonical energy and (v) standard 1+1 energy-class local wellposedness/continuation, all stated explicitly in the theorem documents. The script cannot check infinite-dimensional PDE existence or ensure actual FIN source coefficient binding. No full eigenspectrum digit scan, new external dataset, fitted dynamics, simulator or external hardware was used.

Input and output hash lists are `INPUTS.sha256` and `MANIFEST.sha256`.
