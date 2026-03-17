# N690 Current Strict `T175`: Global Chart Sign Fixing Discharge Theorem (Boundary-Safe, No False-PASS)

Status: `N690_CURRENT_STRICT_T175_GLOBAL_CHART_SIGN_FIXING_DISCHARGE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-17`

## Claim

On the current repo state:

1. the repo exports an explicit chart-level sign-fixing (a `Z2` 0-cochain) rule from strict-core payload weights and the resulting sign-fixed directed representative on `C_v1` (`F690`),
2. applying the same sign-fixing rule to the exported `w_break`-rooted directed representative yields the **same** sign-fixed directed representative (audit `P690`),
3. therefore the sign-fixed directed representative is well-defined (independent of starting exported directed representative) as a tracked `strict_convention` gauge-fixing layer.

This does **not** upgrade directed sign into strict-core physics and does **not** imply kernel-alone/global `QW-2191` discharge.

## Evidence

- `F690` (export of sign-fixed directed representative + explicit sign-fixing data)
- `P690` (independence audit across directed representatives)

## Hard limits

This theorem does **not** claim:

- strict-core physical sign/orientation datum,
- `Aut(Z_12)`-invariant sign canonicity (`N462`),
- kernel-alone/global `QW-2191` discharge,
- operator-level groupoid identities (`N512` boundary),
- ToE closure.

