# Bounded read-only replay

From the repository root, with already installed sympy:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_continuum_phase_lift_research_20261020/verify_research.py --check --output /tmp/fin-continuum-phase-lift-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The 47 controls check scale/hitting algebra, Gaussian radial Jacobians,
the residual integration identity, finite positive constants, all 132
single-copy increment types, actual inherited gamma3, the local dual witness,
and an eight-state variable-clock control. The N144 fixture tests kinematics
only, not MT domain membership, a stationary generator or metastability.
No trajectory simulation, native count enumeration or large LP is run.

`--record` authors only this new package's RESULTS, INPUTS and MANIFEST.
`--check` requires those actual bytes, 13 input hashes and nine predecessor
manifests/166 entries; it writes only the optional /tmp output. No scientific
writer or old finalizer is executed, and current full AGENTS is bound as
the live context rather than equated with a historical full-file hash.

The analytical proof is not mechanically certified by PASS. During drafting,
the variable angle clock was paid by conditional heat-kernel reversal rather
than an unjustified product claim; the quadratic-variation rate was explicitly
stated as mean-square, not L2-norm. The stronger native cosine test then found
a genuine adverse result after the trace test passed. That result is retained;
the diagnostic generator was NOT changed to evade it.

Browser/PDF access is unnecessary for replay. The optional /tmp reading cache
is not a source dependency; its metadata attests actual partial reading,
not proof of native realizability or an original historical download.
