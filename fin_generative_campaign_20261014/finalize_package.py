#!/usr/bin/env python3
"""Attach local provenance, task file hashes, and an archive consistency manifest."""
from pathlib import Path
import os,sys,hashlib,json,shutil
ROOT=Path(__file__).resolve().parent
TOP=ROOT.parent
PLAN=Path('/mnt/data/Wklejony Markdown.md')
assert PLAN.exists(),'source plan not present at mounted supplied path'
shutil.copyfile(PLAN,ROOT/'SOURCE_PLAN_20261014.md')
def hashfile(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def jsonwrite(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
def md(p,s):Path(p).write_text(s.strip()+'\n')
contract=json.loads((ROOT/'RG-C01/MODEL_CONTRACT.json').read_text())
C02=json.loads((ROOT/'RG-C02/RESULTS.json').read_text())
D=json.loads((ROOT/'RG-C03/PROCESS_MAP.json').read_text())
cert=json.loads((ROOT/'RG-C02/INTERVAL_CERTIFICATE.json').read_text())
V=json.loads((ROOT/'RG-C05/PREDICTIONS.json').read_text())
verifier=json.loads((ROOT/'INDEPENDENT_VERIFICATION.json').read_text())
assert verifier['passed']
source_hash=hashfile(ROOT/'SOURCE_PLAN_20261014.md')
md(ROOT/'SOURCE_ACCESS.md',f'''# Source access / Git history limitations

Research date: 2026-10-14. The provided plan is copied unchanged to `SOURCE_PLAN_20261014.md` with local SHA-256 `{source_hash}`. The plan cites baseline **c24773ec**, but that SHA was not verified against GitHub objects here.

## Public remote objects read by web inspection

- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/AGENTS.md — confirms 2026-10-14 guardrails and post-326 reconciliation.
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_post_326_methodology_review/REVIEW_20261014.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_post_326_methodology_review/verification.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-001/DICTIONARY.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-003/PROOF.md
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/COUNTERMODELS.json
- https://raw.githubusercontent.com/hyconiek/Fractal-Nadsoliton-Theory/main/fin_physical_bridge_campaign_01/PHYS-004/CALIBRATION_CONTRACT.json

## Access failures

The GitHub `commits/main` endpoint was restricted in the browsing interface; `git ls-remote` and file downloads in the execution container failed for lack of DNS resolution. Therefore **the latest commit IDs, order, change diff, main HEAD, and exact hash of remote source bytes were not established**. The GitHub homepage listed 1,539 total commits and its root directory was truncated; that does not establish the recent commit sequence. The independently reconstructed A7 implementation faithfully matches the reported 4→3 contrast, but is *not* a byte-verified copy of `safe_seed.py`.

The methodology review and repo guardrails do not show a recoverable FIN-M0 simulator. No Kaggle trajectories or hardware data were imported. This archive is a *new execution package*, not a patched Git checkout. The original repository was not modified.

## Follow-up for a commit-specific audit

From a network-enabled clone, run `git rev-parse HEAD`, `git log -n 20 --date=iso --name-status` and `sha256sum` on the source inputs. Bind the resulting full SHA and source hashes, compare the live source formulas and rerun the producer/verifier. If a formula changed, this handoff cannot be assumed to apply unchanged.
''')
md(ROOT/'README.md',f'''# FIN RG-C01–RG-C05 — executable research handoff

**14 October 2026.** Locally executed bounded finite-N q12 FIN study, based on the supplied plan and public versioned formulas. No repo checkout was available; see `SOURCE_ACCESS.md`.

## What is inside

- `HANDOFF.md`, `CLAIM_REGISTER.json`, `EVIDENCE_LEDGER.md`, `SOURCE_ACCESS.md`.
- `SOURCE_PLAN_20261014.md`: the original user-provided plan, unchanged.
- `RG-C01` to `RG-C05`: reports, structured results, nonconclusions, replay and hashes, plus task-specific artifacts.
- `RG-C02/INTERVAL_CERTIFICATE.json` and `INTERVAL_PROOF.md`: exact-decimal-kernel, finite-alphabet outward enclosure **conditional on mpmath.iv**.
- `RG-C03/PROCESS_MAP.json`, `MEMORY_REPORT.json`: source generator and finite-memory comparisons.
- `RG-C04/COUNTERMODEL_COMPARISON.csv`: normalized countermodel controls.
- `RG-C05/PREREGISTRATION.json`, `PREDICTIONS.json`, `VALIDATION.json`: held-out synthetic model consistency.
- `run_campaign.py`, `interval_certify.py`, `independent_verify.py`: reproducible producer, independent arithmetic certificate, independent second verifier.
- `INDEPENDENT_VERIFICATION.json`: {verifier['number_checks']} separately implemented checks, all passed.
- `RG-C06_CHECKPOINT_FOR_ARCHITECT.md`: decision **not made** by executor; `RG-M00_OPEN_INPUT.md` identifies unavailable inputs.

## Minimal reproduction

Install `numpy`, `scipy`, `mpmath` in a Python 3 environment and run from this directory:

```bash
OPENBLAS_NUM_THREADS=1 python run_campaign.py
OPENBLAS_NUM_THREADS=1 python interval_certify.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
sha256sum -c MANIFEST.sha256
```

**Note:** the first three commands regenerate output files and may alter the checksums in this original manifest because runtime values are recorded; run `sha256sum -c MANIFEST.sha256` **before** replay to verify the distributed archive. After replay regenerate the manifest with the included `finalize_package.py` only after source plan availability is confirmed. `finalize_package.py` requires the original uploaded plan path and is provided for traceability, not needed for running numerical checks.

For manual verification, start with `HANDOFF.md`, `RG-C02/INTERVAL_PROOF.md` and `INDEPENDENT_VERIFICATION.json`.
''')

L=D['FIN_A7'];win='\n'.join('|%.2f|%.9g|%.9g|%.9g|'%(r['t'],r['TV_exact_vs_equilibrium_generator'],r['TV_exact_vs_gibbs_sampler'],r['event00_exact']) for r in L['finite_window'])
md(ROOT/'RG-C03/PROCESS_MAP.md',f'''# Projected CTMC map — FIN_A7

Parent state: 1728 labelled q12 configurations; CTMC heat-bath with rate one per copy. Observation: ordered first two labels (144 states), third hidden. Stationary initial preparation is conditional Gibbs distribution of the hidden label given (0,0); a distinct delta-prepared hidden state tests initial slip. `t` is generator time.

- Source generator row-sum max {L['generator_checks']['row_sum_abs']:.3g}; detailed-balance max {L['generator_checks']['detailed_balance_abs']:.3g}; stationarity max {L['generator_checks']['stationarity_abs']:.3g}.
- Lumpability defect between hidden 0 and hidden 1 at observed (0,0), visible refresh to {L['strong_lumpability_witness']['j_witness']}: {L['strong_lumpability_witness']['maximum_single_visible_refresh_rate_difference']:.9g}. This compares a **real, non-self CTMC jump**.
- Averaged coarse generator B: Eπ[Lf|Y] ; alternative two-label Gibbs sampler D: separate heat bath built from π(Y). The two have the same stationary π(Y) but different transition kernels.

|t|TV(projected, averaged B)|TV(projected, Gibbs sampler D)|P(Y=(0,0)) source|
|---:|---:|---:|---:|
{win}

Reversible equilibrium scalar channel f=1{{Y=(0,0)}} has QLP=v, M0=<v,v>π={L['memory_quadratic_channel']['M0_coupled']:.12g}, M1=<v,Lv>π={L['memory_quadratic_channel']['M1_derivative']:.12g}. Full K(t) in `MEMORY_REPORT.json`. The nonstationary hidden delta-preparation error and its first derivative are **not** contained in a stationary closure; see JSON and `PROCESS_MAP.json`.

These are finite-model numerical statements, not a theory of fundamental physical time.
''')
for task in [f'RG-C0{i}' for i in range(1,6)]:
    path=ROOT/task
    results=json.loads((path/'RESULTS.json').read_text())
    result_hashes={'SOURCE_PLAN_20261014.md':source_hash,'run_campaign.py':hashfile(ROOT/'run_campaign.py')}
    if task!='RG-C01':result_hashes['RG-C01/MODEL_CONTRACT.json']=hashfile(ROOT/'RG-C01/MODEL_CONTRACT.json')
    if task=='RG-C02':result_hashes['interval_certify.py']=hashfile(ROOT/'interval_certify.py')
    if task in ('RG-C04','RG-C05'):
        result_hashes['RG-C02/GENERATED_OPERATORS.json']=hashfile(ROOT/'RG-C02/GENERATED_OPERATORS.json')
        result_hashes['RG-C03/PROCESS_MAP.json']=hashfile(ROOT/'RG-C03/PROCESS_MAP.json')
    results['input_hashes']=result_hashes
    results['source_remote_sha256']='UNAVAILABLE_WEB_ONLY'
    (path/'RESULTS.json').write_text(json.dumps(results,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
    (path/'INPUTS.sha256').write_text(''.join(f'{h}  ../{filename}\n' for filename,h in result_hashes.items()))
    (path/'INPUTS_UNAVAILABLE.md').write_text('# External input limitation\n\nRemote original repository file bytes and latest Git commit objects could not be acquired in the sandbox; therefore exact source-input SHA-256s are OPEN_INPUT. Listed local hashes bind the actual executed reimplementation and provided plan.\n')
# Enhance status and verdict now that interval certificate and standalone second verifier succeeded.
C02p=ROOT/'RG-C02/RESULTS.json';r=json.loads(C02p.read_text());r['arithmetic_status']='binary64 producer; independently outward mpmath.iv exact-decimal source enclosure + rational ANOVA certificate for FIN_A7; secondary verifier passed';r['error_budget']['validated_exact_decimal_kernel_TV_bound_upper']=cert['validated_TV_bound_upper'];r['dependencies']=sorted(set(r['dependencies']+['RG-C02/INTERVAL_CERTIFICATE.json']));jsonwrite(C02p,r)
md(ROOT/'ARITHMETIC_STATUS.md',f'''# Arithmetic and uncertainty boundary

- FIN_A7 static 4→3: all **1728** configurations, exact-decimal kernel `(0.18575,0.1625,1.8)` and gain `5.145228719489142`; outward mpmath.iv input/formula enclosures; exact rational ANOVA on dyadic producer logweights; residual oscillation bounded by interval-projection operator norm. **Validated TV upper:** `{cert['validated_TV_bound_upper']}`. This is conditional on correct `mpmath.iv` implementation and on **this** exact decimal model definition. The source physical A7 remains unproved.
- The true residual ω is inside [`{cert['validated_residual_osc_lower']}`, `{cert['validated_residual_osc_upper']}`]. Full extrema need not be assumed numerically exact; their error is enclosed by 16ε.
- The independent verifier passes {verifier['number_checks']}/{verifier['number_checks']} algebraic/model checks. Independent verifier does **not** constitute formal proof of scipy sparse propagation error.
- Controls' higher-body and TV results are exhaustive **binary64** finite enumerations, not outward-certified against separately defined transcendental control matrices.
- Dynamic CTMC: finite 1728-state sparse generator; reversibility, stationarity and row sums checked numerically. Krylov projected distributions and memory K(t) are numerically propagated without a formal adaptive error envelope.
- RG-C05: static computational holdout uses a producer-based prediction plus a differently enumerated 12^4 parent reference; both depend on the **same** model specification. The rigorous static input uncertainty is covered for FIN only in RG-C02, not an independent data-source error.
- No metrological calibration, device readout, physical time, universal statistical guarantee or external scientific prediction is inferred.
''')
md(ROOT/'REPRODUCTION_STEPS.md','''# Reproduction sequence

1. Inspect `SOURCE_ACCESS.md`; understand that commit provenance is incomplete.
2. Verify package `MANIFEST.sha256` with `sha256sum -c MANIFEST.sha256` in the shipped state (before executing code).
3. Install Python 3 + numpy/scipy/mpmath, run `OPENBLAS_NUM_THREADS=1 python run_campaign.py`.
4. Run `OPENBLAS_NUM_THREADS=1 python interval_certify.py`; do not use an unqualified binary64 bound for exact-kernel formal scope.
5. Run `OPENBLAS_NUM_THREADS=1 python independent_verify.py`, require `passed=true`.
6. Review every `RG-C0x/RESULTS.json`, `NONCONCLUSIONS.md`, and `INPUTS_UNAVAILABLE.md`.
7. If a full repository clone becomes available, bind actual latest Git SHA and original source SHA, inspect diffs, and rerun before attaching findings to upstream.
''')
# report actual result boundaries per branch
md(ROOT/'RG-C04/SPECIFICITY_VERDICT.md',f'''# Countermodel specificity / falsification verdict

All three systems share q12, M4 static reduction, M3 kinetic parent, g={contract['contract']['gain']}, per-copy unit-rate leave-one-out heat-bath and the same 2-label observation. The comparison FULL_POTTS_TRACE and FLAT_P7_TRACE matrices have matching trace to FIN_A7.

|model|third alternating contrast|pair truncation TV|equilibrium scalar memory M0|
|---|---:|---:|---:|
'''+ '\n'.join(f"|{name}|{json.loads((ROOT/'RG-C02/GENERATED_OPERATORS.json').read_text())[2*i]['contrast']:.9f}|{json.loads((ROOT/'RG-C02/GENERATED_OPERATORS.json').read_text())[2*i]['pair_trunc_TV']:.9f}|{D[name]['memory_quadratic_channel']['M0_coupled']:.9f}|" for i,name in enumerate(['FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE']))+'''

**Result:** generated higher-order interactions and failure of strong Markov closure occur in controls too. The existence of these mechanisms therefore **does not establish uniqueness of FIN_A7**. Differences in numerical fingerprints can discriminate *supplied* models, not supply the missing microscopic/physical source for one of them. Further kinetics beyond heat bath are outside this bounded execution.
''')
md(ROOT/'RG-C02/TRUNCATION_CERTIFICATE_SUMMARY.md',f'''# Certified FIN bound summary

- Exact-kernel residual ω ∈ [`{cert['validated_residual_osc_lower']}`, `{cert['validated_residual_osc_upper']}`].
- Hence TV(pair truncated, exact labelled marginal) ≤ `{cert['validated_TV_bound_upper']}` by tanh(ω/4).
- The exact 4→3 formula omits the removed-copy self-diagonal factor because it is an **additive constant in F**; compare raw parent logarithms only up to this constant.
- The theorem is finite and nonspatial. See `INTERVAL_PROOF.md` and `interval_certify.py` for exact-decimal interval arithmetic and dyadic rational proof reduction.
''')
# Root manifest, includes all except itself.
manifest=[]
for p in sorted(ROOT.rglob('*')):
    if not p.is_file() or p.name=='MANIFEST.sha256' or '__pycache__' in p.parts: continue
    manifest.append(f'{hashfile(p)}  {p.relative_to(ROOT).as_posix()}')
(ROOT/'MANIFEST.sha256').write_text('\n'.join(manifest)+'\n')
print('finalized',len(manifest),'files (excluding self)')
