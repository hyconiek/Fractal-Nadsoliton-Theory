# Read-only bounded replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_slow_memory_modulation_20261021/verify_research.py --check --output /tmp/fin-om-20261021.json
git diff --check
```

The new standard-library checker uses exact Fraction/complex-rational
algebra and finite formal Taylor series. It checks shape primitives,
causal convolution identities, slow-time chain-rule factors, exact source
gap margins, initial-slip and resonance controls. Series checks are not
long-time numerical simulations or a mechanical continuum/PDE proof.

NM's existing --check is replayed with isolated /tmp output. Its34 flags
are dependency replay, not additional new OM assertions. No source writer,
new initial packet or external package is run. Old inputs/manifests are
rechecked after execution. --record authors only this new packet;
--check writes only its optional /tmp result. No Internet/PDF cache is
required for the new calculations; NM's optional cache status stays separate.
