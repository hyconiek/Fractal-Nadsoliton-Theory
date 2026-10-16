# Reproducible exact rational certificates: finite symmetric toy networks only

The program `verify_orbit.py` uses Python's `fractions.Fraction` throughout. Its toy models are **not** the actual MT D_N; their role is to detect sign, incidence, capacity, quotient and lift mistakes in the general theorem.

## 36-vertex D6 fixture (six labelled copies of six quotient vertices)

Vertex types: `(level,phase,j)`, levels 0,1,2, phases A,B, cyclic `j in Z/6`. A full reversible native conductance network includes edges from level 0 to level 1: same-phase offsets 0,±1 with weights `1,1/4,1/4`, and cross-phase offsets 0,±1 with weights `1/3,1/6,1/6`. The cross-phase edges from level 0 are removed from the allowed graph H. All four phase-pair edges `(1,p,j)->(2,p',j)` have weight 1 and stay allowed (analogous to a stipulated upper departure layer, **not** a genuine MT boundary).

The exact native demand is `s_(0,p,j)=13/6`, `s_(1,p,j)=-1/6`, `s_(2,p,j)=-2`. It sums to zero. Rotation and reflection of the copy label generate all 12 D6 symmetries.

The exact quotient native-relative bound is `K*=13/9`: the orbit `{(0,A,j):j in Z/6}` has source 13 and allowed outcapacity 9, giving the dual lower bound. A matching construction sends 13 units per phase orbit from level 0 to level 1 along retained arcs, and 6 units from each level-1 phase orbit to **each** level-2 phase orbit. Distribute each quotient arc flow proportionally to its native edge capacities. The lifted 36-vertex flow exactly equals `s`, satisfies `c_e<=13u_e/9`, and has maximum incident conductance `25/6` per state.

With toy reference weights W_x=1 and normalization N=12, the **original vertex-rate cap optimum** is exactly `C*=25/72`. The primal above attains this cap. For a separate exact dual witness take `phi_(l,p,j)=1-l` and `beta_(l,p,j)=1_{l=1}`. Every allowed edge satisfies `phi_tail-phi_head=1<=beta_tail+beta_head`. Then `sum phi*s=50`, `sum beta*W=12`, so `C>=50/(12*12)=25/72`, proving optimality.

Deleting every retained outgoing edge of the level-0 phase-A orbit gives an invariant exact zero-capacity obstruction: this orbit has source 13 but no outgoing allowed edge. This tests the finite `+infinity` branch. This modified graph is an artificial **negative control**, not a legal FIN perturbation.

## Independent complete 6-vertex Z2 cut enumeration

There are 3 levels and two copy indices j=0,1. The full native graph has all cross-copy edges of capacity1 between consecutive levels. The allowed graph removes just the opposite-copy edges at the first layer and keeps all second-layer edges. Its exact source per copy is `(2,0,-2)` across levels. Exhaustive enumeration of all 64 full cuts gives the same sharp `K*=2` as enumeration of just eight orbit-union cuts. This independently checks the reduction on a toy for **all** cuts, not just a few examples.

## Controls and limits

All 96 directed full-network edges and 60 retained edges in the D6 toy are generated independently of the theorem's proof; no floating optimisation is used. The finite examples do not bound actual q12 orbits, do not instantiate the exact A7 Gibbs potential, and do not check the native MT boundary.
