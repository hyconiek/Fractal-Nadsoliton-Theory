# Replay and verification

Run from the package directory:

```bash
python verify.py
sha256sum -c INPUTS.sha256
sha256sum -c MANIFEST.sha256
```

The source hashes in `INPUTS.sha256` are tied to local snapshot file locations relative to `/mnt/data`, so the input check needs the originating mounted conversation files. The manifest covers generated scientific files in this directory and excludes itself and any later ZIP wrapper.

`verify.py` checks B and B inverse exactly in Q(sqrt(3)); exact E2 positivity, spectra, common exit rate and the cubic edge witness; bounded independent matrix-exponential and Poissonization tests; and 120 permutations of the five double rates for the decimal-frozen strict weights. The latter is a finite floating diagnostic; its two surviving permutations are identity and the k1/k5 symmetry relabelling.

Finite matrix controls cannot replace the analytical proof or establish fundamental physics. No simulation of a physical nadsoliton, no new source law, no strong claim about all time-ordered equality patterns.
