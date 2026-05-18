# P1968 Strict P1967 Axes Residual-Z2 Downstream Observable Audit Packet

Status: `P1967_AXES_DOWNSTREAM_Z2_AUDIT_PASS_FOR_PROJECTIVE_OS__SIGN_SENSITIVE_GLOBAL_SELECTOR_STILL_BLOCKED`  
As of: `2026-05-18`

## Goal

Continue after `P1967` by doing the required repo-grep-informed next honest
step:

```text
P1967 gives strict axis-only O(2)->Z2 cuts; now decide which downstream objects
are safe under the remaining residual sign and which still require a directed
sign source.
```

This packet does not reopen the kernel-only obstruction and does not claim full
selector closure.

## Grep result used

The relevant prior work found in `fundamental_action_reconstruction` is:

1. `N493`: residual `Z2` sign flips are gauge/conjugation for the `QW-2190`
   embedding audits.
2. `P709/N706`: Release-7 OS observables actually used downstream are audited
   as residual-sign gauge-irrelevant across all `2^5` pair sign patterns.
3. `P687/N687`: global edge sign coherence is obstructed for directed sections
   under the fixed axis-only transition representatives.
4. `N681/N707`: strict core still lacks a directed/sign-sensitive physical
   orientation datum and the surviving methodology does not discharge global
   `T173` or kernel-alone/global `QW-2191`.

## Symbolic parity lemma

For a residual sign `s in {+1,-1}`, a vector `u`, symmetric quadratic operator
`H`, and linear channel `ell`, the script verifies:

```text
(su)(su)^T = uu^T
(su)^T H (su) = u^T H u
ell^T(su) = s ell^T u
```

Thus:

- projectors/rays and quadratic observables are residual-`Z2` even,
- linear or directed output channels are residual-`Z2` odd unless an additional
  sign convention or strict sign datum is supplied.

## Machine result

The executable probe
`p1968_s918_strict_p1967_axes_residual_z2_downstream_observable_audit.py`
loads the `P1967` axes and existing machine summaries.  It exports a JSON
witness with:

```text
numeric_p1967_projector_and_quadratic_even_pass = true
release7_projective_os_safe_under_residual_z2 = true
global_sign_sensitive_selector_closure_pass = false
global_sign_sensitive_selector_blocked_by_existing_audits = true
```

This is the precise continuation boundary:

```text
projective/quadratic OS continuation is safe;
directed/sign-sensitive global selector closure is still blocked.
```

## No false pass

`P1968` does not claim:

1. directed/sign-sensitive physical orientation datum in strict core,
2. admissible `S_sel_int`,
3. kernel-alone/global `QW-2191` discharge,
4. global `T173` discharge,
5. ToE closure.

## Output artifacts

- Script:
  `p1968_s918_strict_p1967_axes_residual_z2_downstream_observable_audit.py`
- Witness:
  `generated/p1968_s918_strict_p1967_axes_residual_z2_downstream_observable_audit.json`

## Next honest step

If strict ToE closure needs only projective/quadratic OS observables, cite
`P1968` as the residual-`Z2` safety bridge after `P1967`.  If any downstream
observable is genuinely directed/sign-sensitive, build a new global sign
provider or prove a no-go theorem for that directed channel.

## Lay explanation

`P1967` gives axes, but an axis can still be flipped like an arrow with no
arrowhead.  `P1968` checks which quantities notice that flip.  Length-like and
energy-like quantities do not notice it, so the operational OS path is safe.
Arrow-like quantities do notice it, so a full global sign still needs a new
source.
