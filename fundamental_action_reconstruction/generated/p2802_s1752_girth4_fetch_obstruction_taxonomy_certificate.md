# P2802/S1752 girth>=4 fetch obstruction taxonomy certificate

Status: `P2802_GIRTH4_FETCH_OBSTRUCTION_TAXONOMY_CERTIFICATE_NO_IMPORT_NO_CLOSURE`

## Deterministic taxonomy
- candidate_url_count=24
- scheme_counts={'http': 12, 'https': 12}
- https_tunnel_403_count=12
- http_403_count=12
- max_graph6_like_line_count=0
- validated_16828_row_download_count=0
- obstruction_classification=ACCESS_LAYER_BLOCKER_NOT_GRAPH_VALIDATION_FAILURE

## Decision
P2802 classifies P2801 as an access-layer/source-route blocker: all attempted URLs failed before a graph-list body could be validated, so no import or closure is licensed.

## Recommendation
Treat P2802 as an access-layer obstruction certificate.  The next honest step is not another guessed filename probe; it is a working source route that bypasses the 403 layer (manual artifact upload, alternate mirror, source-side access fix, or an installed full generator/toolchain), followed by exact row-count/SHA256/decoding validation and only then quotient/charpoly/complement/orbit auditing.  Otherwise preserve the P2697-P2802 no-canonical-geometry/no-closure certificate.
