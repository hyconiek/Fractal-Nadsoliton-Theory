# FIN 337-WIP — LOW-DEFECT RESPONSE LAW
## Current unfinished research state after task336

Date: 2026-09-28
Status: **WIP — NO PASS/FAIL VERDICT YET**

## 1. Motivation

Task336 stores 2,596 to 11,552 response rows for D<=6 at N=7..10. Under the frozen task333 preparation protocol, however, almost all preparation mass lies much closer to the deep seed.

## 2. Exact D<=2 state count

With eleven non-dominant defect types, the number of weak compositions with total D<=2 is

`1 + 11 + C(12,2) = 78`.

Thus the proposed low-defect basis contains exactly 78 microscopic configurations before any symmetry reduction.

## 3. Mass carried by D<=2 at frozen theta=2

Using the exact task336 state weights and the certified D>6 tail:

| N | D<=6 library size | D<=2 states | certified full mass D<=2 | certified tail D>2 |
|---:|---:|---:|---:|---:|
| 7 | 2,596 | 78 | 99.98806636% | 0.01193364% |
| 8 | 5,648 | 78 | 99.99094605% | 0.00905395% |
| 9 | 9,351 | 78 | 99.99233793% | 0.00766207% |
| 10 | 11,552 | 78 | 99.99307156% | 0.00692844% |

By data processing, replacing the full frozen-protocol preparation by D<=2 alone changes any later observed law by at most the corresponding D>2 tail.

## 4. Proposed response decomposition

At fixed N, write the 13-outcome post-burn response as:

- `R0` for D=0;
- `Delta_a = R(e_a)-R0` for one defect of type a;
- `Delta_ab = R(e_a+e_b)-R(e_a)-R(e_b)+R0` for two-defect interactions (with the obvious repeated-defect convention).

For D<=2 this is an exact algebraic reparameterization of the 78 response rows.

The research question is **not** whether this decomposition exists at fixed N; it does. The question is whether `Delta_a` and `Delta_ab`, after symmetry reduction and mechanically motivated normalization, obey a transferable law in N and later g.

## 5. What has NOT been completed

- no frozen cross-N functional form for Delta_a or Delta_ab;
- no held-out response-law test;
- no cross-g response law;
- no claim that pair interactions are negligible or additive;
- no new-g microscopic opening.

Therefore this must remain WIP.
