# P444 Current Strict T168 (vpsi,g4,g6) Value‑Provider Repo‑Scan Probe

Status: `P444_EXECUTED_CURRENT_STRICT_T168_VPSI_G4_G6_VALUE_PROVIDER_REPO_SCAN_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T168` names the next strict missing provider class beneath the diagonal/local accelerator lane (`T166`):

```text
strict-derived canonical per-site numeric arrays:
  vpsi[0..11], g4[0..11], g6[0..11]
so P437 (N477) can compute the six Sigma_psi{k}_psi{k+6} values and feed P434.
```

`P444` is a repo-state hygiene probe demanded by “no false pass” discipline:

```text
does the repo already contain ANY exported numeric/value object that instantiates
vpsi/g4/g6 in a T168-consumable way (length-12 finite arrays, vpsi nonzero),
and if so, is any such object explicitly marked strict-derived (not toy/probe/extension)?
```

This probe does **not** claim physical admissibility.
It is a search for already-exported numeric instantiations that could make `P437` computable without inventing new
inputs.

## Scan method

Executed by:

- `fundamental_action_reconstruction/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.py`

The scan searches the whole repo (as text) for numeric instantiations of:

1. JSON array-valued assignments:
   - `"vpsi": [<number> ...]`, `"g4": [<number> ...]`, `"g6": [<number> ...]`
2. JSON per-site key assignments (legacy/alt styles):
   - `"vpsi0": <number>`, `"g4_psi0": <number>`, `"g6_psi0": <number>`, etc.
3. Python list/np.array assignments starting with a numeric literal:
   - `vpsi = [<number> ...]`, `g4 = [<number> ...]`, `g6 = [<number> ...]`

Then it attempts to parse candidate JSON files and validates the minimal `T168` computability contract:

1. `vpsi,g4,g6` are present as length-12 lists of finite real numbers,
2. `vpsi[i] != 0` for all `i` (needed by `N477/P437` division premise),
3. candidate is classified only by explicit markers (no silent promotion).

By default, large vendor/cache directories are excluded (same hygiene policy as `P432`).
To scan absolutely everything in a checkout, run with:

```bash
python3 fundamental_action_reconstruction/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.py --full-scan
```

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe.json`
- `fundamental_action_reconstruction/generated/p444_current_strict_t168_vpsi_g4_g6_value_provider_repo_scan_probe_summary.json`

## Result (no false pass)

Update (`2026-03-15`):
after exporting the strict constrained lift/provider `F447` (theorem-anchored by `N483`, provenance-audited by `P448`)
the repo now **does** contain at least one exported `T168`-consumable provider candidate explicitly marked
`strict_derived`.

The scan still finds multiple `T168`-consumable numeric instantiations in `AX` extension-lane artifacts; those remain
useful for computability, but must not be verbally promoted into strict core unless explicitly marked strict-derived
with theorem-level anchors.
