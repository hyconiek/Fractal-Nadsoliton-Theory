# Read-only replay; coefficient arithmetic is not critical-time dynamics

Date: 2026-10-21. From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_two_mode_dissipation_research_20261021/verify_research.py --check
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_subcritical_field_transfer_research_20261021/verify_research.py --check
git diff --check
```

The standard-library Fraction verifier uses exact two-angle Fourier-
Laurent polynomials over Q(i), checking all self/mixed/closed work terms
in two independent amplitude fixtures. It recomputes the rational lower
bounds and their sum-action coercivity, the new Fb profile factor and
normalizations. These do not prove a nonlinear action ODE or a complete
virial inequality. The limiting-response and causal estimates require
their analytical proofs in the two scientific Markdown files.

The own COERCIVITY_CERTIFICATE concerns the FROZEN REFERENCE work form,
not the full nonlinear functional. Old Gamma is reused with its bound
source proof and normalization factor9, not recomputed to more digits.

The script checks bound inputs/protected manifests, own coverage, Python/
JSON syntax and Markdown links/whitespace. It never runs an old source
writer in place. `--record` authors only this fresh packet's metadata.
The optional public PDF/text/screenshot cache in /tmp is not required
for replay and is not copied into the repo. Missing cache is not a
missing scientific source input or reconstructed historical custody.

No browser/network, agent, large simulation, heavy algebra or experiment
is needed for replay. The actual browser reading and precise portions
inspected are recorded separately in LITERATURE_REGISTER.json.
