# P1673 / S623 Helmholtz H4 Boundary + Global Atlas Proxy Packet (strict-only)

Status: `P1673_EXECUTED_H4_BOUNDARY_AND_GLOBAL_ATLAS_PROXY_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel
Domknąć kolejny brak po `P1672/S622`: dodać jawny eksport dla H4
(kontrola członów brzegowych) i spiąć go z globalnym atlasem overlap
w torze odwrotnym `EOM -> L_total`.

## Łańcuch fizyczny
`K_strict -> coeff -> L_SM + L_GR + L_mix -> EOM -> (H1..H4) -> L_total`.

Brak bridge do legacy. Strict-only.

## Eksport
- H4 proxy boundary-balance na trzech klasach warunków brzegowych,
- macierz overlap-consistency dla atlasu `(U0,U1,U2)` jako warunek globalizacji,
- jawna lista niedomkniętych theoremów QG.

## Wynik
Eksport: `generated/p1673_s623_helmholtz_h4_boundary_and_global_atlas_proxy.json`.
Status globalny pozostaje `OPEN_OBLIGATION`.

## Następny uczciwy krok
`S624`: zastąpić proxy H4 formalnym dowodem theorem-level i zintegrować z
pełnym spin-2/spin-0 unitarity theorem + renormalization + background independence.

## Omówienie dla laika
To krok od „lokalnie działa” do „globalnie spójne”. Sprawdzamy, czy teoria
zachowuje się poprawnie także na granicach i przy łączeniu lokalnych kawałków,
ale pełny dowód całości nadal wymaga najtrudniejszych twierdzeń QG.
