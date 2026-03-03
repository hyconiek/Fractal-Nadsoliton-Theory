# External Confirmatory V2 Intake Package

This directory contains templates for a NEW external confirmatory dataset.
Do not reuse design/training datasets.

## Required Files
- `manifest.json` copied from `manifest_template.json` and filled in.
- `pta_v2_pairs.csv` with pair-level PTA features.
- `gw_windows.csv` with GW window-level features.

## Required Protocol Hashes
- PTA V2 protocol hash: `e5bcdc803f5587f790d9c1a70418463ed416760c4fcec72f6cd06b46a92b2f50`
- GW protocol hash: `b9cf21d3d32508e95c6f7cef2e8a953a12f6c7ea7732b6288605c518ab5db5af`

## Minimal Rules
- Data must be external to design dataset used in QW-1848/1850.
- Keep file hashes stable after manifest freeze.
- No post-hoc threshold changes.
