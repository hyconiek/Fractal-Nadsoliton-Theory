# RAPORT QW-1847: PTA REPLICATION TARGET PROPOSAL

- Data UTC: 2026-03-02T22:52:43.733243+00:00
- Strict gate now: PARTIAL / PTA_PROBABILITY_INFERENCE_UNDERPOWERED
- Current n: 14
- Recommended target n: 76
- Additional needed: 62

## Scenario targets for p0=0.90
- Minimal (p_true=1.00): n=29, add=15
- Optimistic (p_true=0.99): n=46, add=32
- Realistic (p_true=0.97): n=76, add=62

## Interpretation
- Przy zachowaniu progu p0=0.90 potrzebna jest istotna eskalacja liczby replikacji, zeby kryterium bylo inferencyjnie domkniete (alpha=0.05, power~0.80).

## Artifacts
- JSON: `report_qw1847_pta_replication_target_proposal.json`
