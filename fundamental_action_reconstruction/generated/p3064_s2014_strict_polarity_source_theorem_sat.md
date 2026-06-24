# P3064/S2014 strict polarity-source theorem SAT certificate

Status: `P3064_STRICT_POLARITY_SOURCE_THEOREM_SAT_NO_CURRENT_EXPORT`

## Finite certificate
- content grep lanes: `3`
- content grep hits: `991`
- primitive atoms: `4`
- SAT rows: `16`
- consistent rows: `9`
- accepted theorem rows: `4`
- minimal accepting atom count: `2`
- minimal accepting rows: `4`
- current exported atoms: `0`
- current row accepted: `False`
- satisfied proof obligations: `3/5`

## Decision
P3064 constructs the missing strict polarity-source theorem as an exhaustive C2/SAT certificate.  A real theorem must contain two primitive non-premise atoms: one absolute source sign and one absolute coupling polarity.  The SAT table enumerates all 16 subsets of the four primitive atoms.  Exactly 4 rows would be accepted if such primitive atoms were exported, and every minimal accepted row has 2 atoms.  The current artifact row has 0 exported primitive atoms and is not accepted, so no G_selector row, QW-2191 discharge, selector closure, L_total, bridge/role transfer, or ToE closure is exported.

## Recommendation
Do not add another partial source/coupling row.  The next proof-grade move must construct exactly one primitive atom from the P3064 SAT theorem: either an absolute non-premise source-sign law or an absolute non-premise coupling-polarity law, then rerun the two-atom acceptance check; otherwise pivot to another P3057 atom with content-first grep.
