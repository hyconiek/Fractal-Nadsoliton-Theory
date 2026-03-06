# RAPORT QW-2107: HZ STRICT DESIGN SEARCH

- Date UTC: 2026-03-04T12:24:47.587050+00:00
- Verdict: **HZ_STRICT_DESIGN_FOUND**
- input_hz_csv: `/home/krzysiek/Pobrane/TOE/edison/external_hz_nodes_qw2099.csv`
- existing_nodes: `3`
- threshold_min_nodes: `5`
- threshold_min_z_span: `0.8`
- threshold_min_e_span: `1.0`
- threshold_max_cond: `8.0`

## Base Metrics
- n_nodes: `3`
- z_span: `0.23`
- e_span: `0.487019`
- design_condition_number: `10.9562`

## Search
- candidate_grid_size: `141`
- k_min_required: `2`
- k_max_searched: `4`
- n_solutions: `4776`

## Top Solutions (z additions)
1. add [0.100, 0.900] -> z_span=0.8, e_span=1.74236, cond=4.39409
2. add [0.100, 0.920] -> z_span=0.82, e_span=1.81136, cond=4.30732
3. add [0.120, 0.920] -> z_span=0.8, e_span=1.78806, cond=4.36056
4. add [0.100, 0.940] -> z_span=0.84, e_span=1.88181, cond=4.22628
5. add [0.120, 0.940] -> z_span=0.82, e_span=1.85852, cond=4.27583
6. add [0.100, 0.960] -> z_span=0.86, e_span=1.95374, cond=4.15094
7. add [0.120, 0.960] -> z_span=0.84, e_span=1.93044, cond=4.19714
8. add [0.140, 0.960] -> z_span=0.82, e_span=1.90629, cond=4.24549
9. add [0.100, 0.980] -> z_span=0.88, e_span=2.02714, cond=4.08122
10. add [0.120, 0.980] -> z_span=0.86, e_span=2.00384, cond=4.12437

## Artifact
- JSON: `report_qw2107_hz_strict_design_search.json`
