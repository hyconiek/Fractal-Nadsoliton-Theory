# N17 Current Existing Kernel Feedback Factorization Obstruction After Current-Pair Chart Reduction Packet Theorem

Status: `N17_DISCHARGED_CURRENT_EXISTING_KERNEL_FEEDBACK_FACTORIZATION_OBSTRUCTION_AFTER_CURRENT_PAIR_CHART_REDUCTION_PACKET_NO_FALSE_PASS`
As of: `2026-03-07`

## Goal

After `R10` and `P14`, the factorization route has become strictly sharper:

- shared provenance is present,
- the explicit current-pair chain and block are present,
- the host carrier is present,
- the host-to-control pushforward is present,
- the chosen current-pair chart reduction is present,
- but the route still does not factorize existing kernel feedback into the
  explicit selector-facing `H3` chain.

`N17` states the strongest honest theorem for that updated route.

## Statement

Consider the route:

```text
existing kernel feedback inside K_total
  + shared frozen-kernel provenance
  + explicit current-pair H3 chain and block
  + host carrier
  + host-to-control pushforward
  + current-pair chart reduction
  -> equivalence/factorization map
```

The theorem is:

> Even after materializing the chosen current-pair chart reduction, the current
> repo still does not identify existing kernel feedback with the explicit
> selector-facing `H3` chain, because the repo still lacks a single
> intertwiner/equality witness identifying the chart-reduced legacy object with
> the computed current-pair `H3` block.

## Result

`N17` discharges:

- a route-specific obstruction theorem for the current repo state after `R10`,
- namely that the route is now reduced to one remaining operator-identification
  witness,
- but still does not factorize existing kernel feedback into the explicit
  selector-facing chain.

## Hard limits

`N17` does not discharge:

- a theorem that no future factorization route can exist,
- a global impossibility theorem for all light-feedback routes,
- `QW-2191`,
- ToE closure.

## Recommended next move

The correct next move is now:

1. materialize the single remaining intertwiner/equality witness,
2. rerun the same factorization route after that addition,
3. keep the theorem negative until that witness is genuinely exported.
