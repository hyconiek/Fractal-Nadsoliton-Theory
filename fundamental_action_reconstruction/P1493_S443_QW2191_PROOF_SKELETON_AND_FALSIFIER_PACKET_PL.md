# P1493 — S4.43 QW-2191 Proof Skeleton + Falsifier Packet (PL)

Status: `P1493_EXECUTED_QW2191_PROOF_SKELETON_AND_FALSIFIER_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Przejść od "kandydata twierdzenia" (`P1492`) do formalnego szkicu dowodu
z jawnie zdefiniowanym falsyfikatorem.

## Decyzja profesorska

Budujemy minimalny szkielet dowodu:

- Lemma L1: bezpieczny zakres `kappa` implikuje `|Delta_SB| <= margin`,
- Lemma L2: w tym zakresie zachodzi `G1 < G0`,
- Lemma L3: orientacja selektora nie zmienia znaku,
- Wniosek C1: istnieje lokalny strict selector-source consistency region.

Równolegle publikujemy falsyfikator:

- dowolny punkt bezpieczny z `G1 >= G0` lub zmianą orientacji
  obala C1.
