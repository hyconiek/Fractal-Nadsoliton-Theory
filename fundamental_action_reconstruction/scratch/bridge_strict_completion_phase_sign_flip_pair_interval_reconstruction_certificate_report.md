# Phase-sign flip-pair interval reconstruction certificate

Status: `ordered-flip-pairs-reconstruct-maximal-node-support-intervals`

## Result

This certificate scans the ordered edge flips from the left anchor and
pairs entry/exit flips to reconstruct the maximal node-support intervals.

## Summary

- Ordered flip edges: `['1->2', '5->6', '7->8', '9->10']`
- Reconstructed intervals: `[{'start': 2, 'end': 5}, {'start': 8, 'end': 9}]`
- Node bits from scan match Z2: `True`
- Edge bits from interval boundaries match Z2: `True`
- Flip count even: `True`
- Flip count equals two components: `True`
- No endpoint support: `True`

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
