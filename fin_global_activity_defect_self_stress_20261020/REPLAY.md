# Bounded read-only replay

From the repository root:

```sh
timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_global_activity_defect_self_stress_20261020/verify_stress.py --check --output /tmp/fin-ad-self-stress-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_global_activity_defect_research_20261020/verify_research.py --check --output /tmp/fin-ad-after-self-stress.json
git diff --check
```

Only standard-library integer/Fraction arithmetic is used. --check reruns
the frozen cases and checks this packet and protected source manifests,
inputs, JSON/AST and local links. --record authors only this new directory.
An explicit /tmp output is the only check-mode write.

A single case may also be inspected with interval_stress.py --case uniform,
interior, face or vertex --output /tmp/<name>.json, under the same 55-second
bound. No new case, gain, N, tolerance or interval subdivision is selected
after an adverse result. Resource exhaustion would be a bounded replay
failure, not a scientific counterexample. These are not stochastic
simulations or independent peer review.
