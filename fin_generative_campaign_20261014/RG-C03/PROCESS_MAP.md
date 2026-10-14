# Projected CTMC map — FIN_A7

Parent state: 1728 labelled q12 configurations; CTMC heat-bath with rate one per copy. Observation: ordered first two labels (144 states), third hidden. Stationary initial preparation is conditional Gibbs distribution of the hidden label given (0,0); a distinct delta-prepared hidden state tests initial slip. `t` is generator time.

- Source generator row-sum max 1.86e-15; detailed-balance max 1.95e-18; stationarity max 2.56e-17.
- Lumpability defect between hidden 0 and hidden 1 at observed (0,0), visible refresh to 1: 0.171914174. This compares a **real, non-self CTMC jump**.
- Averaged coarse generator B: Eπ[Lf|Y] ; alternative two-label Gibbs sampler D: separate heat bath built from π(Y). The two have the same stationary π(Y) but different transition kernels.

|t|TV(projected, averaged B)|TV(projected, Gibbs sampler D)|P(Y=(0,0)) source|
|---:|---:|---:|---:|
|0.15|0.00153232347|0.0200178182|0.950147813|
|0.50|0.0118841589|0.0569494384|0.862192737|
|1.00|0.0304310105|0.0931137884|0.774736539|
|1.50|0.0469030076|0.117588362|0.709148974|

Reversible equilibrium scalar channel f=1{Y=(0,0)} has QLP=v, M0=<v,v>π=0.00995201522562, M1=<v,Lv>π=-0.0180399817805. Full K(t) in `MEMORY_REPORT.json`. The nonstationary hidden delta-preparation error and its first derivative are **not** contained in a stationary closure; see JSON and `PROCESS_MAP.json`.

These are finite-model numerical statements, not a theory of fundamental physical time.
